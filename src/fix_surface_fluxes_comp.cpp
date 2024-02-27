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

#include "fix_surface_fluxes_comp.h"
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

FixSurfaceFluxesComp::FixSurfaceFluxesComp(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR, "Illegal fix surface_fluxes command");
  if (igroup) error->all(FLERR, "fix surface_fluxes must use group all");

  MPI_Comm_rank(world, &me);

  maxbox[0] = 0;
  maxbox[1] = 0;
  maxbox[2] = 0;

  int n = strlen(arg[3]) + 1;
  cordflag = new char[n];
  strcpy(cordflag, arg[3]);

  coordflag = 0;
  binflag = 0;
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
  } else if (strcmp(arg[iarg], "bin") == 0) {
    binflag = 1;
    if (narg != 10) error->all(FLERR, "Illegal fix surface/fluxes command.");
  } else
    error->all(
        FLERR,
        "Illegal fix surface_fluxes command, specify bin, x, y, or z followed by relevant data");

  if (binflag) {
    planecord = atof(arg[6]);
    pinterval = atoi(arg[7]);
    avg_interval = atoi(arg[8]);
    if (me == 0) fp = fopen(arg[9], "w");
  } else {
    planecord = atof(arg[4]);
    pinterval = atoi(arg[5]);
    avg_interval = atoi(arg[6]);
    if (me == 0) fp = fopen(arg[7], "w");
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
      if (binflag)
        sprintf(str, "Cannot open compute surfstressprof file %s", arg[9]);
      else
        sprintf(str, "Cannot open compute surfstressprof file %s", arg[7]);
      error->one(FLERR, str);
    }

    fluxes = NULL;
    surface_sums_2D = NULL;
    surface_sums_2D_global = NULL;

    clearerr(fp);
    fprintf(fp, "# surface-summed data for heat and stress \n");
    if (ferror(fp)) error->one(FLERR, "Error writing file header");
  }
}

/* ---------------------------------------------------------------------- */

FixSurfaceFluxesComp::~FixSurfaceFluxesComp()
{

  if (fp && me == 0) fclose(fp);
  memory->destroy(fluxes);
  memory->destroy(surface_sums_2D_global);
  memory->destroy(surface_sums_2D);
}

/* ---------------------------------------------------------------------- */

int FixSurfaceFluxesComp::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= POST_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesComp::init()
{

  // if non-hybrid, then error if single_enable = 0
  // if hybrid, let hybrid determine if sub-style sets single_enable = 0

  if (pairflag && force->pair == NULL)
    error->all(FLERR, "No pair style defined for fix surface fluxes");
  if (force->pair_match("^hybrid", 0) == NULL && force->pair->single_enable == 0)
    error->all(FLERR, "Pair style does not support fix surface fluxes");
  if (comm->ghost_velocity == 0)
    error->all(FLERR, "fix surface fluxes requires ghost velocity turned on");
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
  // surface_count_x=1;
  //surface_count_y=1;

  memory->create(surface_sums_2D, surface_count_x, surface_count_y, 4, "surfstressprof:2D");
  memory->create(surface_sums_2D_global, surface_count_x, surface_count_y, 4, "surfstressprof:2Dg");
  memory->create(fluxes, 4, "fix_surface_fluxes:fluxes");
  //init to zero
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
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesComp::setup_post_neighbor()
{
  copy_bin_info();
}
/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesComp::post_neighbor()
{
  copy_bin_info();
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesComp::copy_bin_info()
{

  double bbox[3], bsubboxlo[3], bsubboxhi[3];
  double *cutghost = comm->cutghost;
  double *bboxlo = neighbor->bboxlo;
  double *bboxhi = neighbor->bboxhi;
  bsubboxlo[0] = domain->sublo[0] - cutghost[0];
  bsubboxlo[1] = domain->sublo[1] - cutghost[1];
  bsubboxlo[2] = domain->sublo[2] - cutghost[2];
  bsubboxhi[0] = domain->subhi[0] + cutghost[0];
  bsubboxhi[1] = domain->subhi[1] + cutghost[1];
  bsubboxhi[2] = domain->subhi[2] + cutghost[2];
  bbox[0] = bboxhi[0] - bboxlo[0];
  bbox[1] = bboxhi[1] - bboxlo[1];
  bbox[2] = bboxhi[2] - bboxlo[2];
  double cutneighmax = neighbor->cutneighmax;
  double coord;
  //set neighbor bin pointer
  if (neighbor->old_requests[irequest]->copy) {
    int copy_parent = neighbor->old_requests[irequest]->copylist;
    nb = neighbor->neigh_bin[neighbor->old_requests[copy_parent]->index_bin];
  } else
    nb = neighbor->neigh_bin[neighbor->old_requests[irequest]->index_bin];

  //copy bin info
  nbinx = nb->nbinx;
  nbiny = nb->nbiny;
  nbinz = nb->nbinz;
  mbins = nb->mbins;
  mbinx = nb->mbinx;
  mbiny = nb->mbiny;
  mbinz = nb->mbinz;
  mbinxlo = nb->mbinxlo;
  mbinylo = nb->mbinylo;
  mbinzlo = nb->mbinzlo;

  bininvx = nb->bininvx;
  bininvy = nb->bininvy;
  bininvz = nb->bininvz;

  atom2bin = nb->atom2bin;
  bins = nb->bins;
  binhead = nb->binhead;

  //compute high bin limit

  coord = bsubboxhi[0] + SMALL * bbox[0];
  mbinxhi = static_cast<int>((coord - bboxlo[0]) * bininvx);

  coord = bsubboxhi[1] + SMALL * bbox[1];
  mbinyhi = static_cast<int>((coord - bboxlo[1]) * bininvy);

  if (dimension == 3) {
    coord = bsubboxhi[2] + SMALL * bbox[2];
    mbinzhi = static_cast<int>((coord - bboxlo[2]) * bininvz);
  }

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
  sy = static_cast<int>(cutneighmax * bininvy);
  if (sy * binsizey < cutneighmax) sy++;
  sz = static_cast<int>(cutneighmax * bininvz);
  if (sz * binsizez < cutneighmax) sz++;
  if (dimension == 2) sz = 0;

  //determine bins the plane intersects
  int bin_index;
  if (binflag == 0) {
    if (coordflag == 0) {
      bin_index = static_cast<int>((planecord - bboxlo[0]) * bininvx);
      bin_index = MIN(bin_index, nbinx - 1);
      bin_limits[coordflag][0] = bin_index - sx;
      bin_limits[coordflag][1] = bin_index + sx;
      bin_limits[coordflag][1] = MIN(bin_limits[coordflag][1], nbinx - 1);
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
      bin_limits[coordflag][1] = MIN(bin_limits[coordflag][1], nbiny - 1);
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
      bin_limits[coordflag][1] = MIN(bin_limits[coordflag][1], nbinz - 1);
      bin_limits[0][0] = mbinxlo;
      bin_limits[0][1] = mbinxhi;
      bin_limits[1][0] = mbinylo;
      bin_limits[1][1] = mbinyhi;
    }
  }
}
/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesComp::end_of_step()
{
  int k, xx, yy;
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, eng, fpair, factor_coul, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int current_bin, sign;
  int print_flag, average_start, sum_flag;

  double **x = atom->x;
  double **v = atom->v;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  int local_computation = 1;
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
  if (binflag == 0 && sum_flag) {

    for (i = 0; i < atom->nlocal; i++) {
      // skip if atom I is not in either group
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

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx * delx + dely * dely + delz * delz;
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
                (delx * fpair * v[j][0] + dely * fpair * v[j][1] + delz * fpair * v[j][2] +
                 delx * fpair * v[i][0] + dely * fpair * v[i][1] + delz * fpair * v[i][2]) /
                2 * sign;
            surface_sums_2D[0][0][1] -= delx * fpair * sign;
            surface_sums_2D[0][0][2] -= dely * fpair * sign;
            surface_sums_2D[0][0][3] -= delz * fpair * sign;
          } else {
            surface_sums_2D[0][0][0] -=
                (delx * fpair * v[i][0] + dely * fpair * v[i][1] + delz * fpair * v[i][2]) / 2 *
                sign;
            surface_sums_2D[0][0][1] -= delx * fpair * sign / 2;
            surface_sums_2D[0][0][2] -= dely * fpair * sign / 2;
            surface_sums_2D[0][0][3] -= delz * fpair * sign / 2;
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
  //reinit surface accumulation arrays
  if (print_flag) {
    for (i = 0; i < surface_count_x; i++) {
      for (k = 0; k < surface_count_y; k++) {
        for (j = 0; j < 4; j++) {
          surface_sums_2D[i][k][j] = 0.0;
          surface_sums_2D_global[i][k][j] = 0.0;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
