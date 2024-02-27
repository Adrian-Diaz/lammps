/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Adrian Diaz (University of Florida)
------------------------------------------------------------------------- */

#include "fix_surface_plane_fluxes.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "math.h"
#include "memory.h"
#include "modify.h"
#include "nbin.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "nstencil.h"
#include "pair.h"
#include "update.h"
#include "variable.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define DIAGEPSILON 1.0e-6
#define PEPSILON 5.0e-7

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSurfacePlaneFluxes::FixSurfacePlaneFluxes(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  if (narg < 9) error->all(FLERR, "Illegal fix surface_fluxes command");
  if (igroup) error->all(FLERR, "fix surface_fluxes must use group all");

  MPI_Comm_rank(world, &me);

  box_count[0] = 0;
  box_count[1] = 0;
  box_count[2] = 0;

  coordflag = 0;
  tecplot_flag = 0;
  pairflag = 1;
  vector_flag = 1;
  size_vector = 4;
  invoked = -100000000;

  int iarg = 3;
  //specify a single plane or a grid of bins with surfaces
  if (strcmp(arg[iarg], "x") == 0) {
    coordflag = 0;
    if (narg != 9) error->all(FLERR, "Illegal fix surface_fluxes command.");
  } else if (strcmp(arg[iarg], "y") == 0) {
    coordflag = 1;
    if (narg != 9) error->all(FLERR, "Illegal fix surface/fluxes command.");
  } else if (strcmp(arg[iarg], "z") == 0) {
    coordflag = 2;
    if (narg != 9) error->all(FLERR, "Illegal fix surface/fluxes command.");
  } else
    error->all(
        FLERR,
        "Illegal fix surface_fluxes command, specify bin, x, y, or z followed by relevant data");

  planecoord = atof(arg[4]);
  plane_limits[0][0] = atof(arg[5]);
  plane_limits[0][1] = atof(arg[6]);
  plane_limits[1][0] = atof(arg[7]);
  plane_limits[1][1] = atof(arg[8]);

  surface_sums_2D = NULL;
  surface_sums_2D_global = NULL;
  grid_bin_limits = NULL;
  grid_surface_fluxes = NULL;
  global_grid_surface_fluxes = NULL;
  grid_compute_flags = NULL;
}

/* ---------------------------------------------------------------------- */

FixSurfacePlaneFluxes::~FixSurfacePlaneFluxes()
{

  memory->destroy(surface_sums_2D_global);
  memory->destroy(surface_sums_2D);
}

/* ---------------------------------------------------------------------- */

int FixSurfacePlaneFluxes::setmask()
{
  int mask = 0;
  mask |= POST_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSurfacePlaneFluxes::init()
{

  // if non-hybrid, then error if single_enable = 0
  // if hybrid, let hybrid determine if sub-style sets single_enable = 0

  if (pairflag && force->pair == NULL)
    error->all(FLERR, "No pair style defined for fix surface fluxes");
  if (force->pair_match("^hybrid", 0) == NULL && force->pair->single_enable == 0 &&
      force->pair->flux_enable == 0)
    error->all(FLERR, "Pair style does not support fix surface fluxes");
  if (comm->ghost_velocity == 0 && force->newton_pair)
    error->all(FLERR, "fix surface fluxes requires ghost velocity turned on with newton on");
  if (neighbor->dist_check != 1)
    error->all(FLERR, "fix surface fluxes requires neighbor check turned on");
  // error if Kspace style does not compute group/group interactions

  if (pairflag) {
    pair = force->pair;
    cutsq = force->pair->cutsq;
  } else
    pair = NULL;

  dimension = domain->dimension;
  double *prd = domain->prd;
  surface_count_x = 1;
  surface_count_y = 1;
  //store initial timestep
  ptimestep = update->ntimestep;

  //init to zero
  memory->create(surface_sums_2D, surface_count_x, surface_count_y, 4, "fix_surface_fluxes:2D");
  memory->create(surface_sums_2D_global, surface_count_x, surface_count_y, 4,
                 "fix_surface_fluxes:2Dg");
  memory->create(fluxes, 4, "fix_surface_fluxes:fluxes");
  for (int i = 0; i < surface_count_x; i++) {
    for (int k = 0; k < surface_count_y; k++) {
      for (int j = 0; j < 4; j++) {
        surface_sums_2D[i][k][j] = 0.0;
        surface_sums_2D_global[i][k][j] = 0.0;
      }
    }
  }

  irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  if (pair->flux_enable) neighbor->requests[irequest]->newton = 2;

  // bbox lo/hi ptrs = bounding box of entire domain, stored by Domain

  if (domain->triclinic == 0) {
    bboxlo = domain->boxlo;
    bboxhi = domain->boxhi;
  } else {
    bboxlo = domain->boxlo_bound;
    bboxhi = domain->boxhi_bound;
  }

  //test planecoord for feasibility

  if (planecoord > bboxhi[coordflag] || planecoord < bboxlo[coordflag])
    error->all(FLERR,
               "position of the plane for fix surface fluxes must be inside the simulation box");

  debug_sum = 0;
}

/* ---------------------------------------------------------------------- */

void FixSurfacePlaneFluxes::setup_post_neighbor()
{
  bin_info();
}
/* ---------------------------------------------------------------------- */

void FixSurfacePlaneFluxes::post_neighbor()
{
  bin_info();
}

/* ---------------------------------------------------------------------- */

void FixSurfacePlaneFluxes::bin_info()
{

  double bbox[3], bsubboxlo[3], bsubboxhi[3], *sublo, *subhi;
  double grid_point[3], bininv[3];
  int nbin[3], s[3];
  double coord;
  int dim1, dim2, grid_index;
  double *cutghost = comm->cutghost;
  double cutneighmax = neighbor->cutneighmax;
  double skin = neighbor->skin;
  double small = skin * 0.001;
  skin += small;
  sublo = domain->sublo;
  subhi = domain->subhi;
  if (domain->triclinic == 0) {
    bboxlo = domain->boxlo;
    bboxhi = domain->boxhi;
    bsubboxlo[0] = domain->sublo[0] - cutghost[0];
    bsubboxlo[1] = domain->sublo[1] - cutghost[1];
    bsubboxlo[2] = domain->sublo[2] - cutghost[2];
    bsubboxhi[0] = domain->subhi[0] + cutghost[0];
    bsubboxhi[1] = domain->subhi[1] + cutghost[1];
    bsubboxhi[2] = domain->subhi[2] + cutghost[2];
  } else {
    double lo[3], hi[3];
    bboxlo = domain->boxlo_bound;
    bboxhi = domain->boxhi_bound;
    lo[0] = domain->sublo_lamda[0] - cutghost[0];
    lo[1] = domain->sublo_lamda[1] - cutghost[1];
    lo[2] = domain->sublo_lamda[2] - cutghost[2];
    hi[0] = domain->subhi_lamda[0] + cutghost[0];
    hi[1] = domain->subhi_lamda[1] + cutghost[1];
    hi[2] = domain->subhi_lamda[2] + cutghost[2];
    domain->bbox(lo, hi, bsubboxlo, bsubboxhi);
  }
  bbox[0] = bboxhi[0] - bboxlo[0];
  bbox[1] = bboxhi[1] - bboxlo[1];
  bbox[2] = bboxhi[2] - bboxlo[2];

  class NeighRequest *find;
  //set neighbor bin and stencil pointer; first find the parent neigh request if this neigh request is derived

  //find the parent for this pair request if it is constructed as a skip, copy, or halffull list
  find = neighbor->old_requests[irequest];
  int parent_request;
  //see if the neighbor list is derived from another neighbor list (copy,half-full,skip)
  //continue until we find the parent list
  while (find->copy || find->halffull) {
    if (find->copy) parent_request = find->copylist;

    if (find->halffull) parent_request = find->halffulllist;

    find = neighbor->old_requests[parent_request];
  }

  //set neighbor bin and stencil pointers
  nb = neighbor->neigh_bin[find->index_bin];
  ns = neighbor->neigh_stencil[find->index_stencil];

  //copy bin info
  nbin[0] = nbinx = nb->nbinx;
  nbin[1] = nbiny = nb->nbiny;
  nbin[2] = nbinz = nb->nbinz;
  mbins = nb->mbins;
  mbinx = nb->mbinx;
  mbiny = nb->mbiny;
  mbinz = nb->mbinz;
  mbinxlo = nb->mbinxlo;
  mbinylo = nb->mbinylo;
  mbinzlo = nb->mbinzlo;

  bininv[0] = bininvx = nb->bininvx;
  bininv[1] = bininvy = nb->bininvy;
  bininv[2] = bininvz = nb->bininvz;

  atom2bin = nb->atom2bin;
  bins = nb->bins;
  binhead = nb->binhead;

  mbinxhi = mbinx + mbinxlo - 1;
  mbinyhi = mbiny + mbinylo - 1;
  mbinzhi = mbinz + mbinzlo - 1;

  //copy stencil info
  nstencil = ns->nstencil;
  stencil = ns->stencil;
  stencilxyz = ns->stencilxyz;
  nstencil_multi = ns->nstencil_multi_old;
  stencil_multi = ns->stencil_multi_old;
  distsq_multi = ns->distsq_multi_old;
  sx = static_cast<int>(cutneighmax * bininvx);
  if (sx * binsizex < cutneighmax) sx++;
  s[0] = sx;
  sy = static_cast<int>(cutneighmax * bininvy);
  if (sy * binsizey < cutneighmax) sy++;
  s[1] = sy;
  sz = static_cast<int>(cutneighmax * bininvz);
  if (sz * binsizez < cutneighmax) sz++;
  if (dimension == 2) sz = 0;
  s[2] = sz;
  //determine bins the plane intersects
  int bin_index;
  if (coordflag == 0) {
    dim1 = 1;
    dim2 = 2;
    bin_index = static_cast<int>((planecoord - bboxlo[0]) * bininvx);
    bin_index = MIN(bin_index, nbinx - 1);
    bin_limits[coordflag][0] = bin_index - sx;
    bin_limits[coordflag][1] = bin_index + sx;
    //bin_limits[coordflag][1] = MIN(bin_limits[coordflag][1],nbinx-1);
    bin_index = static_cast<int>((plane_limits[0][0] - bboxlo[1]) * bininvy);
    bin_index = MIN(bin_index, nbiny - 1);
    bin_limits[dim1][0] = bin_index - sy;
    bin_index = static_cast<int>((plane_limits[0][1] - bboxlo[1]) * bininvy);
    bin_index = MIN(bin_index, nbiny - 1);
    bin_limits[dim1][1] = bin_index + sy;

    bin_index = static_cast<int>((plane_limits[1][0] - bboxlo[2]) * bininvz);
    bin_index = MIN(bin_index, nbinz - 1);
    bin_limits[dim2][0] = bin_index - sz;
    bin_index = static_cast<int>((plane_limits[1][1] - bboxlo[2]) * bininvz);
    bin_index = MIN(bin_index, nbinz - 1);
    bin_limits[dim2][1] = bin_index + sz;
  }
  if (coordflag == 1) {
    dim1 = 0;
    dim2 = 2;
    bin_index = static_cast<int>((planecoord - bboxlo[1]) * bininvy);
    bin_index = MIN(bin_index, nbiny - 1);
    bin_limits[coordflag][0] = bin_index - sy;
    bin_limits[coordflag][1] = bin_index + sy;
    //bin_limits[coordflag][1] = MIN(bin_limits[coordflag][1],nbiny-1);
    bin_index = static_cast<int>((plane_limits[0][0] - bboxlo[0]) * bininvx);
    bin_index = MIN(bin_index, nbinx - 1);
    bin_limits[dim1][0] = bin_index - sx;
    bin_index = static_cast<int>((plane_limits[0][1] - bboxlo[0]) * bininvx);
    bin_index = MIN(bin_index, nbinx - 1);
    bin_limits[dim1][1] = bin_index + sx;

    bin_index = static_cast<int>((plane_limits[1][0] - bboxlo[2]) * bininvz);
    bin_index = MIN(bin_index, nbinz - 1);
    bin_limits[dim2][0] = bin_index - sz;
    bin_index = static_cast<int>((plane_limits[1][1] - bboxlo[2]) * bininvz);
    bin_index = MIN(bin_index, nbinz - 1);
    bin_limits[dim2][1] = bin_index + sz;
  }
  if (coordflag == 2) {
    dim1 = 0;
    dim2 = 1;
    bin_index = static_cast<int>((planecoord - bboxlo[2]) * bininvz);
    bin_index = MIN(bin_index, nbinz - 1);
    bin_limits[coordflag][0] = bin_index - sz;
    bin_limits[coordflag][1] = bin_index + sz;
    //bin_limits[coordflag][1] = MIN(bin_limits[coordflag][1],nbinz-1);
    bin_index = static_cast<int>((plane_limits[0][0] - bboxlo[0]) * bininvx);
    bin_index = MIN(bin_index, nbinx - 1);
    bin_limits[dim1][0] = bin_index - sx;
    bin_index = static_cast<int>((plane_limits[0][1] - bboxlo[0]) * bininvx);
    bin_index = MIN(bin_index, nbinx - 1);
    bin_limits[dim1][1] = bin_index + sx;

    bin_index = static_cast<int>((plane_limits[1][0] - bboxlo[1]) * bininvy);
    bin_index = MIN(bin_index, nbiny - 1);
    bin_limits[dim2][0] = bin_index - sy;
    bin_index = static_cast<int>((plane_limits[1][1] - bboxlo[1]) * bininvy);
    bin_index = MIN(bin_index, nbiny - 1);
    bin_limits[dim2][1] = bin_index + sy;
  }
  //determine need for local computation of this plane
  local_computation = 1;
  if (bin_limits[0][0] < mbinxlo && bin_limits[0][1] < mbinxlo) local_computation = 0;
  if (bin_limits[1][0] < mbinylo && bin_limits[1][1] < mbinylo) local_computation = 0;
  if (bin_limits[2][0] < mbinzlo && bin_limits[2][1] < mbinzlo) local_computation = 0;
  if (bin_limits[0][1] > mbinxhi && bin_limits[0][0] > mbinxhi) local_computation = 0;
  if (bin_limits[1][1] > mbinyhi && bin_limits[1][0] > mbinyhi) local_computation = 0;
  if (bin_limits[2][1] > mbinzhi && bin_limits[2][0] > mbinzhi) local_computation = 0;

  //loop through bin extent this plane intersects
  if (local_computation) {
    //wrap limits to low and high bin indices of the local task if computation must take place
    bin_limits[0][0] = MAX(bin_limits[0][0], mbinxlo);
    bin_limits[0][1] = MIN(bin_limits[0][1], mbinxhi);
    bin_limits[1][0] = MAX(bin_limits[1][0], mbinylo);
    bin_limits[1][1] = MIN(bin_limits[1][1], mbinyhi);
    bin_limits[2][0] = MAX(bin_limits[2][0], mbinzlo);
    bin_limits[2][1] = MIN(bin_limits[2][1], mbinzhi);
  }
}
/* ---------------------------------------------------------------------- */

double FixSurfacePlaneFluxes::compute_vector(int array_index)
{
  if (invoked != update->ntimestep) {
    int k, xx, yy, zz;
    int i, j, ii, jj, inum, jnum, itype, jtype;
    double xtmp, ytmp, ztmp, delxa[3], proj, intersection_point[3], lparam, delx, dely, delz, *xi,
        *xj, *vi, *vj;
    double vix, viy, viz, vjx, vjy, vjz;
    double rsq, eng, fpair, factor_coul, factor_lj, px, py, pz;
    int *ilist, *jlist, *numneigh, **firstneigh;
    int current_bin, sign, intersection_flag, grid_i[3], grid_j[3], grid_position_i,
        grid_position_j;
    int max_intersection, delta_grid[3], grid_travel[3], current_grid[3], tmp_delta[3], ix[3],
        ix_temp[3];
    int print_flag, average_start, sum_flag;
    double grid_point[3];
    int dim1, dim2, grid_index, grid_index_j, pbc_break;
    double skin = neighbor->skin;
    int *periodicity = domain->periodicity;
    double **x = atom->x;
    double **v = atom->v;
    tagint *molecule = atom->molecule;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    double *special_coul = force->special_coul;
    double *special_lj = force->special_lj;
    int newton_pair = force->newton_pair;
    double interaction_forceij[3], interaction_forceji[3];

    ntimestep = update->ntimestep;
    invoked = 1;
    inum = neighbor->lists[irequest]->inum;
    ilist = neighbor->lists[irequest]->ilist;
    numneigh = neighbor->lists[irequest]->numneigh;
    firstneigh = neighbor->lists[irequest]->firstneigh;

    //initialize sum array for fluxes
    for (i = 0; i < surface_count_x; i++) {
      for (k = 0; k < surface_count_y; k++) {
        for (j = 0; j < 4; j++) {
          surface_sums_2D[i][k][j] = 0.0;
          surface_sums_2D_global[i][k][j] = 0.0;
        }
      }
    }
    //determine if this task needs to accumulate forces to the plane
    if (local_computation) {
      //loop through bin extent this plane intersects
      for (int xbin = bin_limits[0][0]; xbin <= bin_limits[0][1]; xbin++) {
        for (int ybin = bin_limits[1][0]; ybin <= bin_limits[1][1]; ybin++) {
          for (int zbin = bin_limits[2][0]; zbin <= bin_limits[2][1]; zbin++) {
            //compute the current bin id
            current_bin =
                (zbin - mbinzlo) * mbiny * mbinx + (ybin - mbinylo) * mbinx + (xbin - mbinxlo);
            for (i = binhead[current_bin]; i >= 0; i = bins[i]) {
              // skip if atom I is not in either group
              if (i >= atom->nlocal) continue;
              if (!(mask[i] & groupbit)) continue;
              xtmp = x[i][0];
              ytmp = x[i][1];
              ztmp = x[i][2];
              itype = type[i];
              jlist = firstneigh[i];
              jnum = numneigh[i];
              vi = v[i];
              //loop through neighbors of the current atom index from the current bin
              for (jj = 0; jj < jnum; jj++) {
                j = jlist[jj];
                factor_lj = special_lj[sbmask(j)];
                factor_coul = special_coul[sbmask(j)];
                j &= NEIGHMASK;
                xj = x[j];
                vj = v[j];
                // skip if molecule IDs of atoms I,J do not satisfy molflag setting

                delxa[0] = delx = xtmp - x[j][0];
                delxa[1] = dely = ytmp - x[j][1];
                delxa[2] = delz = ztmp - x[j][2];
                rsq = delxa[0] * delxa[0] + delxa[1] * delxa[1] + delxa[2] * delxa[2];
                jtype = type[j];

                if (rsq < cutsq[itype][jtype]) {
                  if (pair->flux_enable) {
                    interaction_forceij[0] = interaction_forceij[1] = interaction_forceij[2] = 0;
                    interaction_forceji[0] = interaction_forceji[1] = interaction_forceji[2] = 0;
                    pair->flux_interaction(i, j, ii, jj, itype, jtype, rsq, interaction_forceij);
                    if (j < atom->nlocal)
                      pair->flux_interaction(j, i, jj, ii, jtype, itype, rsq, interaction_forceji);
                  } else {
                    eng = pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);
                    interaction_forceij[0] = -delx * fpair / 2;
                    interaction_forceij[1] = -dely * fpair / 2;
                    interaction_forceij[2] = -delz * fpair / 2;
                    interaction_forceji[0] = delx * fpair / 2;
                    interaction_forceji[1] = dely * fpair / 2;
                    interaction_forceji[2] = delz * fpair / 2;
                  }
                  if (coordflag == 0) {
                    dim1 = 1;
                    dim2 = 2;
                  }
                  if (coordflag == 1) {
                    dim1 = 0;
                    dim2 = 2;
                  }
                  if (coordflag == 2) {
                    dim1 = 0;
                    dim2 = 1;
                  }

                  //compute perpendicular projection of the line connecting i and j
                  proj = x[i][coordflag] - planecoord;

                  //use the ratio between this projection and the i-j displacement normal to the plane
                  //to define the line parameter (0-1) at the point of intersection
                  lparam = proj / (delxa[coordflag]);
                  if (delxa[coordflag] == 0) continue;

                  //use line parameter to extrapolate the possible intersection point between i-j
                  intersection_point[dim1] = xj[dim1] + delxa[dim1] * (1 - lparam);
                  intersection_point[dim2] = xj[dim2] + delxa[dim2] * (1 - lparam);

                  //test if i-j normal coordinates are on opposing sides of the plane
                  if ((proj < 0 && x[j][coordflag] < planecoord) ||
                      ((proj > 0) && x[j][coordflag] > planecoord))
                    continue;

                  //test the tangential coordinates to determine if the line through i-j crosses the finite sized plane
                  if (intersection_point[dim1] <= plane_limits[0][0] ||
                      intersection_point[dim1] > plane_limits[0][1])
                    continue;
                  if (intersection_point[dim2] <= plane_limits[1][0] ||
                      intersection_point[dim2] > plane_limits[1][1])
                    continue;

                  if (x[i][coordflag] < planecoord) {
                    sign = -1;
                  } else {
                    sign = 1;
                  }
                  surface_sums_2D[0][0][0] +=
                      (interaction_forceij[0] * vi[0] + interaction_forceij[1] * vi[1] +
                       interaction_forceij[2] * vi[2] - interaction_forceji[0] * vj[0] -
                       interaction_forceji[1] * vj[1] - interaction_forceji[2] * vj[2]) *
                      sign;
                  surface_sums_2D[0][0][1] -=
                      (interaction_forceij[0] - interaction_forceji[0]) * sign;
                  surface_sums_2D[0][0][2] -=
                      (interaction_forceij[1] - interaction_forceji[1]) * sign;
                  surface_sums_2D[0][0][3] -=
                      (interaction_forceij[2] - interaction_forceji[2]) * sign;
                }
              }
            }
          }
        }
      }
    }

    MPI_Allreduce(&surface_sums_2D[0][0][0], &surface_sums_2D_global[0][0][0],
                  surface_count_x * surface_count_y * 4, MPI_DOUBLE, MPI_SUM, world);
    invoked = update->ntimestep;
  }
  if (array_index == 0) return surface_sums_2D_global[0][0][0];
  if (array_index == 1) return surface_sums_2D_global[0][0][1];
  if (array_index == 2) return surface_sums_2D_global[0][0][2];
  if (array_index == 3) return surface_sums_2D_global[0][0][3];
}

/* ---------------------------------------------------------------------- */
