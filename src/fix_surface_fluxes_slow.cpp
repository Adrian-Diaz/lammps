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

#include "fix_surface_fluxes_slow.h"
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

#define SMALL 1.0e-6

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSurfaceFluxesSlow::FixSurfaceFluxesSlow(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR, "Illegal fix surface_fluxes command");
  if (igroup) error->all(FLERR, "fix surface_fluxes must use group all");

  MPI_Comm_rank(world, &me);

  box_count[0] = 0;
  box_count[1] = 0;
  box_count[2] = 0;

  int n = strlen(arg[3]) + 1;
  cordflag = new char[n];
  strcpy(cordflag, arg[3]);

  coordflag = 0;
  gridflag = 0;
  tecplot_flag = 0;
  pairflag = 1;

  int iarg = 3;
  //specify a single plane or a grid of bins with surfaces
  if (strcmp(arg[iarg], "x") == 0) {
    coordflag = 0;
    if (narg != 8) error->all(FLERR, "Illegal fix surface_fluxes command.");
  } else if (strcmp(arg[iarg], "y") == 0) {
    coordflag = 1;
    if (narg != 8) error->all(FLERR, "Illegal fix surface/fluxes command.");
  } else if (strcmp(arg[iarg], "z") == 0) {
    coordflag = 2;
    if (narg != 8) error->all(FLERR, "Illegal fix surface/fluxes command.");
  } else if (strcmp(arg[iarg], "grid") == 0) {
    gridflag = 1;
    if (narg < 10 || narg > 11) error->all(FLERR, "Illegal fix surface/fluxes command.");
    if (narg == 11) tecplot_flag = 1;
  } else
    error->all(
        FLERR,
        "Illegal fix surface_fluxes command, specify bin, x, y, or z followed by relevant data");

  if (gridflag) {
    box_count[0] = atoi(arg[4]);
    box_count[1] = atoi(arg[5]);
    box_count[2] = atoi(arg[6]);
    pinterval = atoi(arg[7]);
    avg_interval = atoi(arg[8]);
    if (tecplot_flag && strcmp(arg[9], "tecplot"))
      error->all(FLERR, "Illegal fix surface/fluxes command.");
    if (me == 0) {
      if (tecplot_flag)
        fp = fopen(arg[10], "w");
      else
        fp = fopen(arg[9], "w");
    }
  } else {
    planecord = atof(arg[4]);
    pinterval = atoi(arg[5]);
    avg_interval = atoi(arg[6]);
    if (me == 0) fp = fopen(arg[7], "w");
  }

  //test box counts for feasibility
  if (gridflag) {
    if (box_count[0] <= 0 || box_count[1] <= 0 || box_count[2] <= 0)
      error->all(FLERR, "grid counts for fix surface/fluxes must be greater than zero");
  }

  //test interval settings for feasibility
  if (pinterval <= 0 || avg_interval <= 0)
    error->all(FLERR, "interval settings for fix surface/fluxes must be greater than zero");
  if (pinterval % avg_interval != 0)
    error->all(
        FLERR,
        "print interval must be a multiple of the averaging interval for fix surface/fluxes");
  avg_divisor = pinterval / avg_interval;

  if (me == 0) {
    if (fp == NULL) {
      char str[128];
      if (gridflag)
        sprintf(str, "Cannot open compute surfstressprof file %s", arg[9]);
      else
        sprintf(str, "Cannot open compute surfstressprof file %s", arg[7]);
      error->one(FLERR, str);
    }
    clearerr(fp);
    fprintf(fp, "# surface-summed data for heat and stress \n");
    if (ferror(fp)) error->one(FLERR, "Error writing file header");
  }
  fluxes = NULL;
  surface_sums_2D = NULL;
  surface_sums_2D_global = NULL;
  grid_bin_limits = NULL;
  grid_surface_fluxes = NULL;
  global_grid_surface_fluxes = NULL;
  grid_compute_flags = NULL;
  gather_counts = NULL;
}

/* ---------------------------------------------------------------------- */

FixSurfaceFluxesSlow::~FixSurfaceFluxesSlow()
{

  if (fp && me == 0) fclose(fp);
  memory->destroy(fluxes);
  memory->destroy(surface_sums_2D_global);
  memory->destroy(surface_sums_2D);
  memory->destroy(grid_bin_limits);
  memory->destroy(grid_surface_fluxes);
  memory->destroy(global_grid_surface_fluxes);
  memory->destroy(grid_compute_flags);
}

/* ---------------------------------------------------------------------- */

int FixSurfaceFluxesSlow::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= POST_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesSlow::init()
{

  // if non-hybrid, then error if single_enable = 0
  // if hybrid, let hybrid determine if sub-style sets single_enable = 0

  if (pairflag && force->pair == NULL)
    error->all(FLERR, "No pair style defined for fix surface fluxes");
  if (force->pair_match("^hybrid", 0) == NULL && force->pair->single_enable == 0)
    error->all(FLERR, "Pair style does not support fix surface fluxes");
  if (comm->ghost_velocity == 0 && force->newton_pair)
    error->all(FLERR, "fix surface fluxes requires ghost velocity turned on with newton on");
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
  if (!gridflag) {
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
  }

  irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  // bbox lo/hi ptrs = bounding box of entire domain, stored by Domain

  if (domain->triclinic == 0) {
    bboxlo = domain->boxlo;
    bboxhi = domain->boxhi;
  } else {
    bboxlo = domain->boxlo_bound;
    bboxhi = domain->boxhi_bound;
  }

  //test planecord for feasibility
  if (!gridflag) {
    if (planecord > bboxhi[coordflag] || planecord < bboxlo[coordflag])
      error->all(FLERR,
                 "position of the plane for fix surface fluxes must be inside the simulation box");
  }

  //compute grid box sizes
  if (gridflag) {
    double bbox[3];
    bbox[0] = bboxhi[0] - bboxlo[0];
    bbox[1] = bboxhi[1] - bboxlo[1];
    bbox[2] = bboxhi[2] - bboxlo[2];
    grid_size[0] = bbox[0] / box_count[0];
    grid_size[1] = bbox[1] / box_count[1];
    grid_size[2] = bbox[2] / box_count[2];
  }

  //allocate grid arrays
  //init grid to zero
  if (gridflag) {
    memory->create(grid_bin_limits, box_count[0] * box_count[1] * box_count[2], 6, 3, 2,
                   "fix_surface_fluxes:grid_bin_limits");
    memory->create(grid_surface_fluxes, box_count[0] * box_count[1] * box_count[2], 3, 4,
                   "fix_surface_fluxes:grid_bin_limits");
    memory->create(global_grid_surface_fluxes, box_count[0] * box_count[1] * box_count[2], 3, 4,
                   "fix_surface_fluxes:grid_bin_limits");
    memory->create(grid_compute_flags, box_count[0] * box_count[1] * box_count[2], 3,
                   "fix_surface_fluxes:grid_compute_flags");
    for (int ix = 0; ix < box_count[0]; ix++)
      for (int iy = 0; iy < box_count[1]; iy++)
        for (int iz = 0; iz < box_count[2]; iz++)
          for (int is = 0; is < 3; is++)
            for (int ifl = 0; ifl < 4; ifl++) {
              grid_surface_fluxes[box_count[0] * box_count[1] * iz + box_count[0] * iy + ix][is]
                                 [ifl] = 0;
              global_grid_surface_fluxes[box_count[0] * box_count[1] * iz + box_count[0] * iy + ix]
                                        [is][ifl] = 0;
            }
  }

  debug_sum = 0;
  //allocate gather counts array for local grid data
  memory->create(gather_counts, comm->nprocs, "fix_surface_fluxes:gather_counts");
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesSlow::setup_post_neighbor()
{
  bin_info();
}
/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesSlow::post_neighbor()
{
  bin_info();
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesSlow::bin_info()
{

  double bbox[3], bsubboxlo[3], bsubboxhi[3], *sublo, *subhi;
  double grid_point[3], bininv[3];
  int nbin[3], s[3];
  double coord;
  double plane_limits[2][2];
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

  //set neighbor bin pointer
  if (neighbor->old_requests[irequest]->copy) {
    int copy_parent = neighbor->old_requests[irequest]->copylist;
    nb = neighbor->neigh_bin[neighbor->old_requests[copy_parent]->index_bin];
  } else
    nb = neighbor->neigh_bin[neighbor->old_requests[irequest]->index_bin];

  if (gridflag) {
    grid_size[0] = bbox[0] / box_count[0];
    grid_size[1] = bbox[1] / box_count[1];
    grid_size[2] = bbox[2] / box_count[2];
  }
  /*
  //check if grid size is smaller than bin size; currently unsupported
  if(gridflag){
  if (grid_size[0]<=nb->binsizex)
    error->one(FLERR,"fix surface fluxes requires the grid size to be larger than the bin size; which is usually cutoff/2");
    if (grid_size[1]<=nb->binsizey)
    error->one(FLERR,"fix surface fluxes requires the grid size to be larger than the bin size, which is usually cutoff/2");
    if (grid_size[2]<=nb->binsizez)
    error->one(FLERR,"fix surface fluxes requires the grid size to be larger than the bin size; which is usually cutoff/2");
  }
  */
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

  //set neighbor stencil pointer
  if (neighbor->old_requests[irequest]->copy) {
    int copy_parent = neighbor->old_requests[irequest]->copylist;
    ns = neighbor->neigh_stencil[neighbor->old_requests[copy_parent]->index_stencil];
  } else
    ns = neighbor->neigh_stencil[neighbor->old_requests[irequest]->index_stencil];

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
  if (gridflag == 0) {
    if (coordflag == 0) {
      bin_index = static_cast<int>((planecord - bboxlo[0]) * bininvx);
      bin_index = MIN(bin_index, nbinx - 1);
      bin_limits[coordflag][0] = bin_index - sx;
      bin_limits[coordflag][1] = bin_index + sx;
      //bin_limits[coordflag][1] = MIN(bin_limits[coordflag][1],nbinx-1);
      bin_limits[1][0] = mbinylo;
      bin_limits[1][1] = mbinyhi;
      bin_limits[2][0] = mbinzlo;
      bin_limits[2][1] = mbinzhi;
    }
    if (coordflag == 1) {
      bin_index = static_cast<int>((planecord - bboxlo[1]) * bininvy);
      bin_index = MIN(bin_index, nbiny - 1);
      bin_limits[coordflag][0] = bin_index - sy;
      bin_limits[coordflag][1] = bin_index + sy;
      //bin_limits[coordflag][1] = MIN(bin_limits[coordflag][1],nbiny-1);
      bin_limits[0][0] = mbinxlo;
      bin_limits[0][1] = mbinxhi;
      bin_limits[2][0] = mbinzlo;
      bin_limits[2][1] = mbinzhi;
    }
    if (coordflag == 2) {
      bin_index = static_cast<int>((planecord - bboxlo[2]) * bininvz);
      bin_index = MIN(bin_index, nbinz - 1);
      bin_limits[coordflag][0] = bin_index - sz;
      bin_limits[coordflag][1] = bin_index + sz;
      //bin_limits[coordflag][1] = MIN(bin_limits[coordflag][1],nbinz-1);
      bin_limits[0][0] = mbinxlo;
      bin_limits[0][1] = mbinxhi;
      bin_limits[1][0] = mbinylo;
      bin_limits[1][1] = mbinyhi;
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
  } else {
    //determine need for local flux accumulation of each grid surface
    //compute grid index limits for this subbox
    for (int idim = 0; idim < 3; idim++) {
      sub_grid_limits[idim][0] =
          static_cast<int>((bsubboxlo[idim] - bboxlo[idim]) / grid_size[idim]);
      sub_grid_limits[idim][0] = MAX(0, sub_grid_limits[idim][0]);
      sub_grid_limits[idim][1] =
          static_cast<int>((bsubboxhi[idim] - bboxlo[idim]) / grid_size[idim]);
      sub_grid_limits[idim][1] = MIN(box_count[idim] - 1, sub_grid_limits[idim][1]);
    }
    //loop over grid surfaces to determine bin extents of each
    for (int ix = sub_grid_limits[0][0]; ix <= sub_grid_limits[0][1]; ix++)
      for (int iy = sub_grid_limits[1][0]; iy <= sub_grid_limits[1][1]; iy++)
        for (int iz = sub_grid_limits[2][0]; iz <= sub_grid_limits[2][1]; iz++)
          for (int is = 0; is < 3; is++) {
            //compute current grid surface coordination and size
            if (is == 0) {
              dim1 = 1;
              dim2 = 2;
            }
            if (is == 1) {
              dim1 = 0;
              dim2 = 2;
            }
            if (is == 2) {
              dim1 = 0;
              dim2 = 1;
            }
            grid_index = box_count[0] * box_count[1] * iz + box_count[0] * iy + ix;
            grid_point[0] = grid_size[0] * ix + bboxlo[0];
            grid_point[1] = grid_size[1] * iy + bboxlo[1];
            grid_point[2] = grid_size[2] * iz + bboxlo[2];
            planecord = grid_point[is] + grid_size[is];
            plane_limits[0][0] = grid_point[dim1];
            plane_limits[0][1] = grid_point[dim1] + grid_size[dim1];
            plane_limits[1][0] = grid_point[dim2];
            plane_limits[1][1] = grid_point[dim2] + grid_size[dim2];
            //augment plane limits by skin/2 if this plane borders global simulation box (since particles can leave the box by at most skin/2)
            if (is == 0) {
              if (iy == 0) plane_limits[0][0] -= skin / 2;
              if (iz == 0) plane_limits[1][0] -= skin / 2;
              if (iy == box_count[1] - 1) plane_limits[0][1] += skin / 2;
              if (iz == box_count[2] - 1) plane_limits[1][1] += skin / 2;
            }
            if (is == 1) {
              if (ix == 0) plane_limits[0][0] -= skin / 2;
              if (iz == 0) plane_limits[1][0] -= skin / 2;
              if (ix == box_count[0] - 1) plane_limits[0][1] += skin / 2;
              if (iz == box_count[2] - 1) plane_limits[1][1] += skin / 2;
            }
            if (is == 2) {
              if (ix == 0) plane_limits[0][0] -= skin / 2;
              if (iy == 0) plane_limits[1][0] -= skin / 2;
              if (ix == box_count[0] - 1) plane_limits[0][1] += skin / 2;
              if (iy == box_count[1] - 1) plane_limits[1][1] += skin / 2;
            }

            if (planecord >= bboxhi[is])
              bin_index = static_cast<int>((planecord - bboxhi[is]) * bininv[is]) + nbin[is];
            else if (planecord >= bboxlo[is]) {
              bin_index = static_cast<int>((planecord - bboxlo[is]) * bininv[is]);
              bin_index = MIN(bin_index, nbin[is] - 1);
            } else
              bin_index = static_cast<int>((planecord - bboxlo[is]) * bininv[is]) - 1;

            grid_bin_limits[grid_index][is][is][0] = bin_index - s[is];
            grid_bin_limits[grid_index][is][is][1] = bin_index + s[is];

            if (plane_limits[0][0] >= bboxhi[dim1])
              bin_index =
                  static_cast<int>((plane_limits[0][0] - bboxhi[dim1]) * bininv[dim1]) + nbin[dim1];
            else if (plane_limits[0][0] >= bboxlo[dim1]) {
              bin_index = static_cast<int>((plane_limits[0][0] - bboxlo[dim1]) * bininv[dim1]);
              bin_index = MIN(bin_index, nbin[dim1] - 1);
            } else
              bin_index = static_cast<int>((plane_limits[0][0] - bboxlo[dim1]) * bininv[dim1]) - 1;

            grid_bin_limits[grid_index][is][dim1][0] = bin_index - s[dim1];

            if (plane_limits[0][1] >= bboxhi[dim1])
              bin_index =
                  static_cast<int>((plane_limits[0][1] - bboxhi[dim1]) * bininv[dim1]) + nbin[dim1];
            else if (plane_limits[0][1] >= bboxlo[dim1]) {
              bin_index = static_cast<int>((plane_limits[0][1] - bboxlo[dim1]) * bininv[dim1]);
              bin_index = MIN(bin_index, nbin[dim1] - 1);
            } else
              bin_index = static_cast<int>((plane_limits[0][1] - bboxlo[dim1]) * bininv[dim1]) - 1;

            grid_bin_limits[grid_index][is][dim1][1] = bin_index + s[dim1];

            if (plane_limits[1][0] >= bboxhi[dim2])
              bin_index =
                  static_cast<int>((plane_limits[1][0] - bboxhi[dim2]) * bininv[dim2]) + nbin[dim2];
            else if (plane_limits[1][0] >= bboxlo[dim2]) {
              bin_index = static_cast<int>((plane_limits[1][0] - bboxlo[dim2]) * bininv[dim2]);
              bin_index = MIN(bin_index, nbin[dim2] - 1);
            } else
              bin_index = static_cast<int>((plane_limits[1][0] - bboxlo[dim2]) * bininv[dim2]) - 1;

            grid_bin_limits[grid_index][is][dim2][0] = bin_index - s[dim2];

            if (plane_limits[1][1] >= bboxhi[dim2])
              bin_index =
                  static_cast<int>((plane_limits[1][1] - bboxhi[dim2]) * bininv[dim2]) + nbin[dim2];
            else if (plane_limits[1][1] >= bboxlo[dim2]) {
              bin_index = static_cast<int>((plane_limits[1][1] - bboxlo[dim2]) * bininv[dim2]);
              bin_index = MIN(bin_index, nbin[dim2] - 1);
            } else
              bin_index = static_cast<int>((plane_limits[1][1] - bboxlo[dim2]) * bininv[dim2]) - 1;

            grid_bin_limits[grid_index][is][dim2][1] = bin_index + s[dim2];
          }

    //loop over subbox grid limits to determine if grid surfaces require local computation
    for (int ix = sub_grid_limits[0][0]; ix <= sub_grid_limits[0][1]; ix++)
      for (int iy = sub_grid_limits[1][0]; iy <= sub_grid_limits[1][1]; iy++)
        for (int iz = sub_grid_limits[2][0]; iz <= sub_grid_limits[2][1]; iz++)
          for (int is = 0; is < 3; is++) {
            grid_index = box_count[0] * box_count[1] * iz + box_count[0] * iy + ix;
            //determine need for local computation of this plane
            local_computation = 1;
            if (grid_bin_limits[grid_index][is][0][0] < mbinxlo &&
                grid_bin_limits[grid_index][is][0][1] < mbinxlo)
              local_computation = 0;
            if (grid_bin_limits[grid_index][is][1][0] < mbinylo &&
                grid_bin_limits[grid_index][is][1][1] < mbinylo)
              local_computation = 0;
            if (grid_bin_limits[grid_index][is][2][0] < mbinzlo &&
                grid_bin_limits[grid_index][is][2][1] < mbinzlo)
              local_computation = 0;
            if (grid_bin_limits[grid_index][is][0][1] > mbinxhi &&
                grid_bin_limits[grid_index][is][0][0] > mbinxhi)
              local_computation = 0;
            if (grid_bin_limits[grid_index][is][1][1] > mbinyhi &&
                grid_bin_limits[grid_index][is][1][0] > mbinyhi)
              local_computation = 0;
            if (grid_bin_limits[grid_index][is][2][1] > mbinzhi &&
                grid_bin_limits[grid_index][is][2][0] > mbinzhi)
              local_computation = 0;
            grid_compute_flags[grid_index][is] = local_computation;
            //loop through bin extent this plane intersects
            if (local_computation) {
              //wrap limits to low and high bin indices of the local task if computation must take place
              grid_bin_limits[grid_index][is][0][0] =
                  MAX(grid_bin_limits[grid_index][is][0][0], mbinxlo);
              grid_bin_limits[grid_index][is][0][1] =
                  MIN(grid_bin_limits[grid_index][is][0][1], mbinxhi);
              grid_bin_limits[grid_index][is][1][0] =
                  MAX(grid_bin_limits[grid_index][is][1][0], mbinylo);
              grid_bin_limits[grid_index][is][1][1] =
                  MIN(grid_bin_limits[grid_index][is][1][1], mbinyhi);
              grid_bin_limits[grid_index][is][2][0] =
                  MAX(grid_bin_limits[grid_index][is][2][0], mbinzlo);
              grid_bin_limits[grid_index][is][2][1] =
                  MIN(grid_bin_limits[grid_index][is][2][1], mbinzhi);
            }
          }
    //define grid counts for grid data
    int grid_size = (sub_grid_limits[0][1] - sub_grid_limits[0][0] + 1) *
        (sub_grid_limits[1][1] - sub_grid_limits[1][0] + 1) *
        (sub_grid_limits[2][1] - sub_grid_limits[2][0] + 1);
    gather_counts[comm->me] = grid_size * 3 * 4;
  }
}
/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesSlow::end_of_step()
{
  int k, xx, yy, zz;
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx[3], proj, intersection_point[3], lparam;
  double rsq, eng, fpair, factor_coul, factor_lj, px, py, pz;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int current_bin, sign, intersection_flag;
  int print_flag, average_start, sum_flag;
  double grid_point[3];
  double plane_limits[2][2];
  int dim1, dim2, grid_index;
  double skin = neighbor->skin;

  double **x = atom->x;
  double **v = atom->v;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  ntimestep = update->ntimestep;

  inum = pair->list->inum;
  ilist = pair->list->ilist;
  numneigh = pair->list->numneigh;
  firstneigh = pair->list->firstneigh;
  print_flag = 0;
  //determine whether to print on this timestep
  if (ntimestep % pinterval == 0) print_flag = 1;
  if (sum_flag % avg_interval == 0) sum_flag = 1;
  //check if simulation started on a timestep that is not a multiple of the print interval
  if ((ptimestep) % pinterval != 0) {
    //perform no computation until the timestep succeeds a multiple of the print interval
    if ((ntimestep - ptimestep - 1) / pinterval < 1) {
      print_flag = 0;
      sum_flag = 0;
    }
  }
  if (gridflag == 0 && sum_flag) {

    /*
        if(pair->manybody_flag==1){
            //pair->single_three(1,1,1,1,1,1,1,fluxes,planecords,coordflag);
        }
        if(pair->manybody_flag==1){
          surface_sums_2D[0][0][0]+=fluxes[0];
          surface_sums_2D[0][0][1]+=fluxes[1];
          surface_sums_2D[0][0][2]+=fluxes[2];
          surface_sums_2D[0][0][3]+=fluxes[3];
        }
        */

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
              //loop through neighbors of the current atom index from the current bin
              for (jj = 0; jj < jnum; jj++) {
                j = jlist[jj];
                factor_lj = special_lj[sbmask(j)];
                factor_coul = special_coul[sbmask(j)];
                j &= NEIGHMASK;

                // skip if molecule IDs of atoms I,J do not satisfy molflag setting

                delx[0] = xtmp - x[j][0];
                delx[1] = ytmp - x[j][1];
                delx[2] = ztmp - x[j][2];
                rsq = delx[0] * delx[0] + delx[1] * delx[1] + delx[2] * delx[2];
                jtype = type[j];

                if (rsq < cutsq[itype][jtype]) {
                  eng = pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

                  // energy only computed once so tally full amount
                  // force tally is jgroup acting on igroup

                  if ((x[i][coordflag] > planecord && x[j][coordflag] > planecord) ||
                      (x[i][coordflag] < planecord && x[j][coordflag] < planecord))
                    continue;    // skip if pair trajectory
                  //does not intersect plane

                  if (x[i][coordflag] < planecord) {
                    sign = -1;
                  } else {
                    sign = 1;
                  }
                  if (newton_pair) {
                    surface_sums_2D[0][0][0] -=
                        (delx[0] * fpair * v[j][0] + delx[1] * fpair * v[j][1] +
                         delx[2] * fpair * v[j][2] + delx[0] * fpair * v[i][0] +
                         delx[1] * fpair * v[i][1] + delx[2] * fpair * v[i][2]) /
                        2 * sign;
                    surface_sums_2D[0][0][1] -= delx[0] * fpair * sign;
                    surface_sums_2D[0][0][2] -= delx[1] * fpair * sign;
                    surface_sums_2D[0][0][3] -= delx[2] * fpair * sign;
                  } else {
                    surface_sums_2D[0][0][0] -=
                        (delx[0] * fpair * v[i][0] + delx[1] * fpair * v[i][1] +
                         delx[2] * fpair * v[i][2]) /
                        2 * sign;
                    surface_sums_2D[0][0][1] -= delx[0] * fpair * sign / 2;
                    surface_sums_2D[0][0][2] -= delx[1] * fpair * sign / 2;
                    surface_sums_2D[0][0][3] -= delx[2] * fpair * sign / 2;
                  }
                }
              }
            }
          }
        }
      }
    }
    if (print_flag)
      MPI_Allreduce(&surface_sums_2D[0][0][0], &surface_sums_2D_global[0][0][0],
                    surface_count_x * surface_count_y * 4, MPI_DOUBLE, MPI_SUM, world);

    if (me == 0) {
      if (print_flag) {
        clearerr(fp);
        fprintf(fp, "#Timestep  Number-of-Grids-X  Number-of-Grids-Y \n");

        fprintf(fp, " %d %d %d\n", ntimestep, surface_count_x, surface_count_y);

        for (xx = 0; xx < surface_count_x; xx++) {
          for (yy = 0; yy < surface_count_y; yy++) {
            fprintf(fp, "%d %d %g %g %g %g \n", xx + 1, yy + 1,
                    surface_sums_2D_global[xx][yy][0] / avg_divisor,
                    surface_sums_2D_global[xx][yy][1] / avg_divisor,
                    surface_sums_2D_global[xx][yy][2] / avg_divisor,
                    surface_sums_2D_global[xx][yy][3] / avg_divisor);
          }
        }

        fflush(fp);
      }
    }
  }

  //grid flux calculation
  if (gridflag == 1 && sum_flag) {
    //loop through bin extent this plane intersects
    for (int ix = sub_grid_limits[0][0]; ix <= sub_grid_limits[0][1]; ix++)
      for (int iy = sub_grid_limits[1][0]; iy <= sub_grid_limits[1][1]; iy++)
        for (int iz = sub_grid_limits[2][0]; iz <= sub_grid_limits[2][1]; iz++)
          for (int is = 0; is < 3; is++) {
            coordflag = is;
            grid_index = box_count[0] * box_count[1] * iz + box_count[0] * iy + ix;
            if (grid_compute_flags[grid_index][is]) {
              if (is == 0) {
                dim1 = 1;
                dim2 = 2;
              }
              if (is == 1) {
                dim1 = 0;
                dim2 = 2;
              }
              if (is == 2) {
                dim1 = 0;
                dim2 = 1;
              }
              grid_point[0] = grid_size[0] * ix + bboxlo[0];
              grid_point[1] = grid_size[1] * iy + bboxlo[1];
              grid_point[2] = grid_size[2] * iz + bboxlo[2];
              planecord = grid_point[is] + grid_size[is];
              plane_limits[0][0] = grid_point[dim1];
              plane_limits[0][1] = grid_point[dim1] + grid_size[dim1];
              plane_limits[1][0] = grid_point[dim2];
              plane_limits[1][1] = grid_point[dim2] + grid_size[dim2];
              //augment plane limits by skin/2 if this plane borders global simulation box (since particles can leave the box by at most skin/2)
              if (is == 0) {
                if (iy == 0) plane_limits[0][0] -= skin / 2;
                if (iz == 0) plane_limits[1][0] -= skin / 2;
                if (iy == box_count[1] - 1) plane_limits[0][1] += skin / 2;
                if (iz == box_count[2] - 1) plane_limits[1][1] += skin / 2;
              }
              if (is == 1) {
                if (ix == 0) plane_limits[0][0] -= skin / 2;
                if (iz == 0) plane_limits[1][0] -= skin / 2;
                if (ix == box_count[0] - 1) plane_limits[0][1] += skin / 2;
                if (iz == box_count[2] - 1) plane_limits[1][1] += skin / 2;
              }
              if (is == 2) {
                if (ix == 0) plane_limits[0][0] -= skin / 2;
                if (iy == 0) plane_limits[1][0] -= skin / 2;
                if (ix == box_count[0] - 1) plane_limits[0][1] += skin / 2;
                if (iy == box_count[1] - 1) plane_limits[1][1] += skin / 2;
              }

              for (int xbin = grid_bin_limits[grid_index][is][0][0];
                   xbin <= grid_bin_limits[grid_index][is][0][1]; xbin++) {
                for (int ybin = grid_bin_limits[grid_index][is][1][0];
                     ybin <= grid_bin_limits[grid_index][is][1][1]; ybin++) {
                  for (int zbin = grid_bin_limits[grid_index][is][2][0];
                       zbin <= grid_bin_limits[grid_index][is][2][1]; zbin++) {
                    //compute the current bin id
                    current_bin = (zbin - mbinzlo) * mbiny * mbinx + (ybin - mbinylo) * mbinx +
                        (xbin - mbinxlo);
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
                      //loop through neighbors of the current atom index from the current bin
                      for (jj = 0; jj < jnum; jj++) {
                        j = jlist[jj];
                        factor_lj = special_lj[sbmask(j)];
                        factor_coul = special_coul[sbmask(j)];
                        j &= NEIGHMASK;

                        delx[0] = xtmp - x[j][0];
                        delx[1] = ytmp - x[j][1];
                        delx[2] = ztmp - x[j][2];
                        rsq = delx[0] * delx[0] + delx[1] * delx[1] + delx[2] * delx[2];
                        jtype = type[j];
                        intersection_flag = 1;
                        //compute perpendicular projection of the line connecting i and j
                        proj = x[i][coordflag] - planecord;

                        //use the ratio between this projection and the i-j displacement normal to the plane
                        //to define the line parameter (0-1) at the point of intersection
                        lparam = proj / (delx[coordflag]);
                        if (delx[coordflag] == 0) intersection_flag = 0;

                        //use line parameter to extrapolate the possible intersection point between i-j
                        intersection_point[dim1] = x[j][dim1] + delx[dim1] * (1 - lparam);
                        intersection_point[dim2] = x[j][dim2] + delx[dim2] * (1 - lparam);

                        //test if i-j normal coordinates are on opposing sides of the plane
                        if ((proj < 0 && x[j][coordflag] < planecord) ||
                            ((proj > 0) && x[j][coordflag] > planecord))
                          intersection_flag = 0;
                        if (intersection_flag != 0 &&
                            (intersection_point[dim1] <= plane_limits[0][0] ||
                             intersection_point[dim1] > plane_limits[0][1])) {
                          int outcount = 1;
                        }
                        if (intersection_flag != 0 &&
                            (intersection_point[dim2] <= plane_limits[1][0] ||
                             intersection_point[dim2] > plane_limits[1][1])) {
                          int outcount = 1;
                        }
                        //test the tangential coordinates to determine if the line through i-j crosses the finite sized plane
                        if (intersection_point[dim1] <= plane_limits[0][0] ||
                            intersection_point[dim1] > plane_limits[0][1])
                          intersection_flag = 0;
                        if (intersection_point[dim2] <= plane_limits[1][0] ||
                            intersection_point[dim2] > plane_limits[1][1])
                          intersection_flag = 0;

                        if (!intersection_flag) continue;    // skip if pair trajectory
                        //does not intersect plane

                        if (rsq < cutsq[itype][jtype]) {
                          eng =
                              pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

                          // energy only computed once so tally full amount
                          // force tally is jgroup acting on igroup

                          if (x[i][coordflag] < planecord) {
                            sign = -1;
                          } else {
                            sign = 1;
                          }
                          if (newton_pair) {
                            grid_surface_fluxes[grid_index][is][0] -=
                                (delx[0] * fpair * v[j][0] + delx[1] * fpair * v[j][1] +
                                 delx[2] * fpair * v[j][2] + delx[0] * fpair * v[i][0] +
                                 delx[1] * fpair * v[i][1] + delx[2] * fpair * v[i][2]) /
                                2 * sign;
                            grid_surface_fluxes[grid_index][is][1] -= delx[0] * fpair * sign;
                            grid_surface_fluxes[grid_index][is][2] -= delx[1] * fpair * sign;
                            grid_surface_fluxes[grid_index][is][3] -= delx[2] * fpair * sign;
                          } else {
                            grid_surface_fluxes[grid_index][is][0] -=
                                (delx[0] * fpair * v[i][0] + delx[1] * fpair * v[i][1] +
                                 delx[2] * fpair * v[i][2]) /
                                2 * sign;
                            grid_surface_fluxes[grid_index][is][1] -= delx[0] * fpair * sign / 2;
                            grid_surface_fluxes[grid_index][is][2] -= delx[1] * fpair * sign / 2;
                            grid_surface_fluxes[grid_index][is][3] -= delx[2] * fpair * sign / 2;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }

    if (print_flag)
      MPI_Reduce(&grid_surface_fluxes[0][0][0], &global_grid_surface_fluxes[0][0][0],
                 box_count[0] * box_count[1] * box_count[2] * 4 * 3, MPI_DOUBLE, MPI_SUM, 0, world);
    if (me == 0) {
      if (print_flag) {
        if (!tecplot_flag) {
          clearerr(fp);
          fprintf(fp, "#Timestep  Grid Count-X  Grid Count-Y Grid Count-Z \n");

          fprintf(fp, " %d %d %d %d\n", ntimestep, box_count[0], box_count[1], box_count[2]);

          for (xx = 0; xx < box_count[0]; xx++) {
            for (yy = 0; yy < box_count[1]; yy++) {
              for (zz = 0; zz < box_count[2]; zz++) {
                for (int is = 0; is < domain->dimension; is++) {
                  grid_index = box_count[0] * box_count[1] * zz + box_count[0] * yy + xx;
                  fprintf(fp, "%d %d %d %d %g %g %g %g \n", xx + 1, yy + 1, zz + 1, is + 1,
                          global_grid_surface_fluxes[grid_index][is][0] / avg_divisor,
                          global_grid_surface_fluxes[grid_index][is][1] / avg_divisor,
                          global_grid_surface_fluxes[grid_index][is][2] / avg_divisor,
                          global_grid_surface_fluxes[grid_index][is][3] / avg_divisor);
                }
              }
            }
          }
          //debug_sum
          for (int sumy = 0; sumy < box_count[1]; sumy++)
            for (int sumz = 0; sumz < box_count[2]; sumz++) {
              grid_index = box_count[0] * box_count[1] * sumz + box_count[0] * sumy + 0;
              debug_sum += global_grid_surface_fluxes[grid_index][0][1] / avg_divisor;
            }
          fprintf(fp, " \n \n \n");
          fprintf(fp, "DEBUG SUM  %g \n", debug_sum);
          fprintf(fp, " \n \n \n");
        } else {
          clearerr(fp);
          fprintf(fp,
                  "TITLE=\"Surface Flux Results\""
                  "\n");
          fprintf(fp,
                  "VARIABLES = \"x\", \"y\", \"z\", \"q\", \"sx\", \"sy\", \"sz\""
                  "\n");
          fprintf(fp,
                  "ZONE T=\"load step %d\", I=%d, J=2, DATAPACKING=POINT"
                  "\n",
                  ntimestep, box_count[0] * box_count[1] * box_count[2] * 3 / 2);

          for (xx = 0; xx < box_count[0]; xx++) {
            for (yy = 0; yy < box_count[1]; yy++) {
              for (zz = 0; zz < box_count[2]; zz++) {
                grid_index = box_count[0] * box_count[1] * zz + box_count[0] * yy + xx;
                grid_point[0] = grid_size[0] * xx + bboxlo[0];
                grid_point[1] = grid_size[1] * yy + bboxlo[1];
                grid_point[2] = grid_size[2] * zz + bboxlo[2];
                for (int is = 0; is < domain->dimension; is++) {
                  if (is == 0) {
                    dim1 = 1;
                    dim2 = 2;
                    px = grid_point[is] + grid_size[is];
                    py = grid_point[dim1] + grid_size[dim1] / 2;
                    pz = grid_point[dim2] + grid_size[dim2] / 2;
                  }
                  if (is == 1) {
                    dim1 = 0;
                    dim2 = 2;
                    py = grid_point[is] + grid_size[is];
                    px = grid_point[dim1] + grid_size[dim1] / 2;
                    pz = grid_point[dim2] + grid_size[dim2] / 2;
                  }
                  if (is == 2) {
                    dim1 = 0;
                    dim2 = 1;
                    pz = grid_point[is] + grid_size[is];
                    px = grid_point[dim1] + grid_size[dim1] / 2;
                    py = grid_point[dim2] + grid_size[dim2] / 2;
                  }
                  fprintf(fp, "%g %g %g %g %g %g %g \n", px, py, pz,
                          global_grid_surface_fluxes[grid_index][is][0] / avg_divisor,
                          global_grid_surface_fluxes[grid_index][is][1] / avg_divisor,
                          global_grid_surface_fluxes[grid_index][is][2] / avg_divisor,
                          global_grid_surface_fluxes[grid_index][is][3] / avg_divisor);
                }
              }
            }
          }
          //debug_sum
          for (int sumy = 0; sumy < box_count[1]; sumy++)
            for (int sumz = 0; sumz < box_count[2]; sumz++) {
              grid_index = box_count[0] * box_count[1] * sumz + box_count[0] * sumy + 0;
              debug_sum += global_grid_surface_fluxes[grid_index][0][1] / avg_divisor;
            }
          fprintf(fp, " \n \n \n");
          fprintf(fp, "DEBUG SUM  %g \n", debug_sum);
          fprintf(fp, " \n \n \n");
        }
        fflush(fp);
      }
    }
  }

  //reinit surface accumulation arrays
  if (print_flag && !gridflag) {
    for (i = 0; i < surface_count_x; i++) {
      for (k = 0; k < surface_count_y; k++) {
        for (j = 0; j < 4; j++) {
          surface_sums_2D[i][k][j] = 0.0;
          surface_sums_2D_global[i][k][j] = 0.0;
        }
      }
    }
  }
  if (print_flag && gridflag) {
    debug_sum = 0;
    for (int ix = sub_grid_limits[0][0]; ix <= sub_grid_limits[0][1]; ix++)
      for (int iy = sub_grid_limits[1][0]; iy <= sub_grid_limits[1][1]; iy++)
        for (int iz = sub_grid_limits[2][0]; iz <= sub_grid_limits[2][1]; iz++)
          for (int is = 0; is < 3; is++)
            for (int ifl = 0; ifl < 4; ifl++) {
              grid_surface_fluxes[box_count[0] * box_count[1] * iz + box_count[0] * iy + ix][is]
                                 [ifl] = 0;
              global_grid_surface_fluxes[box_count[0] * box_count[1] * iz + box_count[0] * iy + ix]
                                        [is][ifl] = 0;
            }
  }
}

/* ---------------------------------------------------------------------- */
