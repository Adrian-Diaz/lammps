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

#include "fix_surface_fluxes_atom_op.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
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

#define DIAGEPSILON 1.0e-7
#define PEPSILON 2.0e-7
#define DELEPSILON 2.0e-7
#define MAXBINCONTENT 100    //used to bound local memory
#define EXPAND 100
#define MAXBINOVERLAP 10

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSurfaceFluxesAtomOp::FixSurfaceFluxesAtomOp(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  if (narg != 6 && narg != 9) error->all(FLERR, "Illegal fix surface_fluxes command");

  MPI_Comm_rank(world, &me);

  peratom_flag = 1;
  size_peratom_cols = 24;
  comm_reverse = 24;
  peratom_freq = 1;
  dynamic_group_allow = 1;

  box_size[0] = 0;
  box_size[1] = 0;
  box_size[2] = 0;
  box_center[0] = 0;
  box_center[1] = 0;
  box_center[2] = 0;

  pairflag = 1;

  box_size[0] = atof(arg[3]);
  box_size[1] = atof(arg[4]);
  box_size[2] = atof(arg[5]);

  if (narg > 6) {
    box_center[0] = atof(arg[6]);
    box_center[1] = atof(arg[7]);
    box_center[2] = atof(arg[8]);
  }

  //test box sizes for feasibility
  if (box_size[0] <= 0 || box_size[1] <= 0 || box_size[2] <= 0)
    error->all(FLERR, "box sizes for fix surface/fluxes/atom must be greater than zero");

  //make sure box_center displacement is smaller than half the box size (so that the atom stays inside the box)
  if ((box_center[0] > 0 && box_size[0] <= box_center[0]) ||
      (box_center[0] < 0 && -box_size[0] >= box_center[0]))
    error->all(FLERR,
               "box center cannot be larger than or equal to half the box size (atom must remain "
               "inside the box)");
  if ((box_center[1] > 0 && box_size[1] <= box_center[1]) ||
      (box_center[1] < 0 && -box_size[1] >= box_center[1]))
    error->all(FLERR,
               "box center cannot be larger than or equal to half the box size (atom must remain "
               "inside the box)");
  if ((box_center[2] > 0 && box_size[2] <= box_center[2]) ||
      (box_center[2] < 0 && -box_size[2] >= box_center[2]))
    error->all(FLERR,
               "box center cannot be larger than or equal to half the box size (atom must remain "
               "inside the box)");

  memory->create(bin_overlap_limits, 6, "nbin_CAC:bin_overlap_limits");
  maxbin = 0;
  nmax = 0;
  max_intersect_count = 0;
  fluxes = NULL;
  group_include = NULL;
  atom_fluxes = NULL;
  bin_ncontent = NULL;
  bin_content = NULL;
  bin_expansion_counts = NULL;
  ij_plane_intersects = NULL;
  scan_flags = NULL;
  first_alloc = 0;
  max_bin_expansion_count = 0;
}

/* ---------------------------------------------------------------------- */

FixSurfaceFluxesAtomOp::~FixSurfaceFluxesAtomOp()
{

  memory->destroy(fluxes);
  memory->destroy(atom_fluxes);
  memory->destroy(bin_ncontent);
  memory->destroy(bin_expansion_counts);
  memory->destroy(group_include);

  for (int i = 0; i < maxbin; i++) memory->destroy(bin_content[i]);
  memory->sfree(bin_content);
  memory->destroy(ij_plane_intersects);
}

/* ---------------------------------------------------------------------- */

int FixSurfaceFluxesAtomOp::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  mask |= POST_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesAtomOp::init()
{

  // if non-hybrid, then error if single_enable = 0
  // if hybrid, let hybrid determine if sub-style sets single_enable = 0

  if (pairflag && force->pair == NULL)
    error->all(FLERR, "No pair style defined for fix surface fluxes");
  if (force->pair_match("^hybrid", 0) == NULL && force->pair->single_enable == 0 &&
      force->pair->flux_enable == 0)
    error->all(FLERR, "Pair style does not support fix surface fluxes");
  if (comm->ghost_velocity == 0)
    error->all(FLERR, "fix surface fluxes requires ghost velocity turned on");
  if (neighbor->dist_check != 1)
    error->warning(
        FLERR, "fix surface fluxes requires might not be accurate with neighbor check turned off");
  // error if Kspace style does not compute group/group interactions

  if (pairflag) {
    pair = force->pair;
    cutsq = force->pair->cutsq;
  } else
    pair = NULL;

  dimension = domain->dimension;
  double *prd = domain->prd;

  //init to zero
  array_atom =
      memory->grow(atom_fluxes, atom->nlocal + atom->nghost, 24, "fix_surface_fluxes:atom_fluxes");
  for (int i = 0; i < atom->nlocal + atom->nghost; i++) {
    for (int k = 0; k < 24; k++) { atom_fluxes[i][k] = 0; }
  }

  int newton_pair = force->newton_pair;
  /*
  irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  if(!newton_pair){
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  }
  */
  // bbox lo/hi ptrs = bounding box of entire domain, stored by Domain

  if (domain->triclinic == 0) {
    bboxlo = domain->boxlo;
    bboxhi = domain->boxhi;
  } else {
    bboxlo = domain->boxlo_bound;
    bboxhi = domain->boxhi_bound;
  }

  if (pairflag) {
    irequest = neighbor->request(this, instance_me);
    neighbor->requests[irequest]->pair = 0;
    neighbor->requests[irequest]->fix = 1;
    if (pair->flux_enable) neighbor->requests[irequest]->newton = 2;
  }
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesAtomOp::setup_post_neighbor()
{
  array_atom =
      memory->grow(atom_fluxes, atom->nlocal + atom->nghost, 24, "fix_surface_fluxes:atom_fluxes");
  memory->grow(scan_flags, atom->nlocal + atom->nghost, "scan_flags");
  //init scan flags
  for (int init = 0; init < atom->nlocal + atom->nghost; init++) scan_flags[init] = 0;
  bin_info();
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesAtomOp::setup(int vflag)
{
  compute_flux();
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesAtomOp::post_neighbor()
{
  array_atom =
      memory->grow(atom_fluxes, atom->nlocal + atom->nghost, 24, "fix_surface_fluxes:atom_fluxes");
  memory->grow(scan_flags, atom->nlocal + atom->nghost, "scan_flags");
  for (int init = 0; init < atom->nlocal + atom->nghost; init++) scan_flags[init] = 0;
  bin_info();
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesAtomOp::bin_info()
{
  double skin = neighbor->skin;
  double *cutghost = comm->cutghost;
  double cutneighmax = neighbor->cutneighmax;

  //check if box sizes exceed twice the cutoff radius; this would exceed a communication limit
  if (box_size[0] / 2 > cutghost[0] - skin || box_size[1] / 2 > cutghost[1] - skin ||
      box_size[2] / 2 > cutghost[2] - skin)
    error->all(FLERR,
               "per atom box sizes for fix surface/fluxes/atom must be smaller than twice the "
               "cutoff radius");

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

  binsize[0] = binsizex = nb->binsizex;
  binsize[1] = binsizey = nb->binsizey;
  binsize[2] = binsizez = nb->binsizez;

  atom2bin = nb->atom2bin;
  bins = nb->bins;
  binhead = nb->binhead;

  mbinxhi = mbinx + mbinxlo - 1;
  mbinyhi = mbiny + mbinylo - 1;
  mbinzhi = mbinz + mbinzlo - 1;

  bboxlo = neighbor->bboxlo;
  bboxhi = neighbor->bboxhi;

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
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesAtomOp::end_of_step()
{
  compute_flux();
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesAtomOp::compute_flux()
{
  double xtmp, ytmp, ztmp, delxa[3], proj, intersection_point[3], lparam, delx, dely, delz, *xi,
      *xj, intersect_count;
  double vix, viy, viz, vjx, vjy, vjz;
  double rsq, eng, fpair, factor_coul, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int sign, intersection_flag, intersection_flag2, coordflag, normal_flag;
  int *periodicity = domain->periodicity;
  double **x = atom->x;
  double **v = atom->v;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double interaction_forceij[3], interaction_forceji[3];
  double grid_point[3];
  double bounding_boxlo[3], bounding_boxhi[3];
  double skin = neighbor->skin;
  int current_bin, intersection_flags[3];
  int i, j, ii, jj, inum, jnum, itype, jtype, i_bin, j_bin, bi, mindir;
  double plane_limits[2][2], planecoord, delta, test_max, dels[3], smin;
  int dim1, dim2, grid_index;
  double *cutghost = comm->cutghost;
  double cutneighmax = neighbor->cutneighmax;
  int max_intersection, delta_grid[3], grid_travel[3], current_grid[3], tmp_delta[3], ix[3],
      ix_temp[3];
  int ibin[3], jbin[3];

  inum = neighbor->lists[irequest]->inum;
  ilist = neighbor->lists[irequest]->ilist;
  numneigh = neighbor->lists[irequest]->numneigh;
  firstneigh = neighbor->lists[irequest]->firstneigh;

  //initialize fluxes array
  for (i = 0; i < atom->nlocal + atom->nghost; i++) {
    for (int k = 0; k < 24; k++) { atom_fluxes[i][k] = 0; }
  }

  //determine and store bin overlaps for each atom's box
  //start by allocating bin arrays for atom boxes

  if (mbins > maxbin) allocate_arrays(mbins, 0);

  //initialize bin counts
  for (int init = 0; init < maxbin; init++) bin_ncontent[init] = 0;

  //loop through atom boxes and decide which bins they overlap with
  for (i = 0; i < nlocal + atom->nghost; i++) {
    box2bins(i);
    for (int overlapx = bin_overlap_limits[0]; overlapx <= bin_overlap_limits[3]; overlapx++) {
      for (int overlapy = bin_overlap_limits[1]; overlapy <= bin_overlap_limits[4]; overlapy++) {
        for (int overlapz = bin_overlap_limits[2]; overlapz <= bin_overlap_limits[5]; overlapz++) {
          current_bin = (overlapz - mbinzlo) * mbiny * mbinx + (overlapy - mbinylo) * mbinx +
              (overlapx - mbinxlo);
          //check bin memory allocation is large enough and grow if needed
          if (bin_ncontent[current_bin] ==
              MAXBINCONTENT + bin_expansion_counts[current_bin] * EXPAND) {
            bin_expansion_counts[current_bin]++;
            memory->grow(bin_content[current_bin],
                         MAXBINCONTENT + bin_expansion_counts[current_bin] * EXPAND,
                         "neigh:bin_content grow");
          }
          bin_content[current_bin][bin_ncontent[current_bin]++] = i;
        }
      }
    }
  }

  inum = neighbor->lists[irequest]->inum;
  ilist = neighbor->lists[irequest]->ilist;
  numneigh = neighbor->lists[irequest]->numneigh;
  firstneigh = neighbor->lists[irequest]->firstneigh;

  //loop through neighbors of each atom to determine the set of surfaces intersecting with each pair's connecting line
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xi = x[i];
    i_bin = atom2bin[i];

    if (x[i][0] >= bboxhi[0])
      ibin[0] = static_cast<int>((x[i][0] - bboxhi[0]) * bininvx) + nbinx;
    else if (x[i][0] >= bboxlo[0]) {
      ibin[0] = static_cast<int>((x[i][0] - bboxlo[0]) * bininvx);
      ibin[0] = min(ibin[0], nbinx - 1);
    } else
      ibin[0] = static_cast<int>((x[i][0] - bboxlo[0]) * bininvx) - 1;

    if (x[i][1] >= bboxhi[1])
      ibin[1] = static_cast<int>((x[i][1] - bboxhi[1]) * bininvy) + nbiny;
    else if (x[i][1] >= bboxlo[1]) {
      ibin[1] = static_cast<int>((x[i][1] - bboxlo[1]) * bininvy);
      ibin[1] = min(ibin[1], nbiny - 1);
    } else
      ibin[1] = static_cast<int>((x[i][1] - bboxlo[1]) * bininvy) - 1;

    if (x[i][2] >= bboxhi[2])
      ibin[2] = static_cast<int>((x[i][2] - bboxhi[2]) * bininvz) + nbinz;
    else if (x[i][2] >= bboxlo[2]) {
      ibin[2] = static_cast<int>((x[i][2] - bboxlo[2]) * bininvz);
      ibin[2] = min(ibin[2], nbinz - 1);
    } else
      ibin[2] = static_cast<int>((x[i][2] - bboxlo[2]) * bininvz) - 1;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    vix = v[i][0];
    viy = v[i][1];
    viz = v[i][2];

    for (jj = 0; jj < jnum; jj++) {
      ij_nplane_intersects = 0;
      j = jlist[jj];
      j &= NEIGHMASK;
      xj = x[j];
      j_bin = atom2bin[j];
      /*loop from the starting bin of the line to the terminating bin and search each bins contents
      for intersecting planes of other atom boxes */

      delxa[0] = delx = xtmp - x[j][0];
      delxa[1] = dely = ytmp - x[j][1];
      delxa[2] = delz = ztmp - x[j][2];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      vjx = v[j][0];
      vjy = v[j][1];
      vjz = v[j][2];
      rsq = delxa[0] * delxa[0] + delxa[1] * delxa[1] + delxa[2] * delxa[2];
      jtype = type[j];

      //loop through the set of possible intersecting atom boxes
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

        if (interaction_forceij[0] == 0 && interaction_forceij[1] == 0 &&
            interaction_forceij[2] == 0 && interaction_forceji[0] == 0 &&
            interaction_forceji[1] == 0 && interaction_forceji[2] == 0)
          continue;

        grid_travel[0] = 0;
        grid_travel[1] = 0;
        grid_travel[2] = 0;

        if (x[j][0] >= bboxhi[0])
          jbin[0] = static_cast<int>((x[j][0] - bboxhi[0]) * bininvx) + nbinx;
        else if (x[j][0] >= bboxlo[0]) {
          jbin[0] = static_cast<int>((x[j][0] - bboxlo[0]) * bininvx);
          jbin[0] = min(jbin[0], nbinx - 1);
        } else
          jbin[0] = static_cast<int>((x[j][0] - bboxlo[0]) * bininvx) - 1;

        if (x[j][1] >= bboxhi[1])
          jbin[1] = static_cast<int>((x[j][1] - bboxhi[1]) * bininvy) + nbiny;
        else if (x[j][1] >= bboxlo[1]) {
          jbin[1] = static_cast<int>((x[j][1] - bboxlo[1]) * bininvy);
          jbin[1] = min(jbin[1], nbiny - 1);
        } else
          jbin[1] = static_cast<int>((x[j][1] - bboxlo[1]) * bininvy) - 1;

        if (x[j][2] >= bboxhi[2])
          jbin[2] = static_cast<int>((x[j][2] - bboxhi[2]) * bininvz) + nbinz;
        else if (x[j][2] >= bboxlo[2]) {
          jbin[2] = static_cast<int>((x[j][2] - bboxlo[2]) * bininvz);
          jbin[2] = min(jbin[2], nbinz - 1);
        } else
          jbin[2] = static_cast<int>((x[j][2] - bboxlo[2]) * bininvz) - 1;

        delta_grid[0] = ibin[0] - jbin[0];
        delta_grid[1] = ibin[1] - jbin[1];
        delta_grid[2] = ibin[2] - jbin[2];

        /*check if delta_grid is due to numerical error in positons;
      //can lead to spurious calculations of intersection points later
      //that break the bin grid search since delx is small. These differences negligibly effect
      the calculation.
      */
        for (int delta_check = 0; delta_check < domain->dimension; delta_check++) {
          if (ibin[delta_check] != jbin[delta_check]) {
            if ((x[i][delta_check] - x[j][delta_check] > 0) &&
                (x[i][delta_check] - x[j][delta_check] < DELEPSILON))
              delta_grid[delta_check]--;
            if ((x[i][delta_check] - x[j][delta_check] < 0) &&
                (x[i][delta_check] - x[j][delta_check] > -DELEPSILON))
              delta_grid[delta_check]++;
          }
        }

        tmp_delta[0] = delta_grid[0];
        tmp_delta[1] = delta_grid[1];
        tmp_delta[2] = delta_grid[2];
        //check signs
        if (delta_grid[0] < 0) tmp_delta[0] = -delta_grid[0];
        if (delta_grid[1] < 0) tmp_delta[1] = -delta_grid[1];
        if (delta_grid[2] < 0) tmp_delta[2] = -delta_grid[2];
        //max number of intersections determined by sum of cartesian displacements in grid
        max_intersection = tmp_delta[0] + tmp_delta[1] + tmp_delta[2];

        //search bins the connecting line passes through for atom box intersections
        for (int iloop = 0; iloop < max_intersection + 1; iloop++) {
          //determine the next bin being intersected as we traverse the line from i to j
          //starts at bin occupied by xi
          ix[0] = ibin[0] + grid_travel[0];
          ix[1] = ibin[1] + grid_travel[1];
          ix[2] = ibin[2] + grid_travel[2];
          grid_index =
              (ix[2] - mbinzlo) * mbiny * mbinx + (ix[1] - mbinylo) * mbinx + (ix[0] - mbinxlo);

          //search the current bin for possible atom box intersections
          for (int bin_count = 0; bin_count < bin_ncontent[grid_index]; bin_count++) {
            bi = bin_content[grid_index][bin_count];
            if (scan_flags[bi]) continue;
            if ((bi == i || bi == j) && (mask[bi] & groupbit)) {
              //check if max_intersect for this pair has been reached and reallocate accordingly
              if (ij_nplane_intersects == max_intersect_count) {
                max_intersect_count += MAXBINOVERLAP;
                memory->grow(ij_plane_intersects, max_intersect_count,
                             "FixSurfaceFluxesAtom:ij_plane_intersects");
              }
              ij_plane_intersects[ij_nplane_intersects++] = bi;
              scan_flags[bi] = 1;
              intersect_count = 0;
              //determine which of the 6 planes of the atom box are intersected by this i-j pair
              for (int isl = 0; isl < 2 * domain->dimension; isl++) {
                int is = isl / 2;
                if (isl % 2 == 0)
                  normal_flag = 1;
                else
                  normal_flag = -1;
                coordflag = is;
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

                //test negative and positive sides of the box dimension
                if (isl % 2 == 0)
                  planecoord = x[bi][is] - box_size[is] / 2 + box_center[is];
                else
                  planecoord = x[bi][is] + box_size[is] / 2 + box_center[is];
                plane_limits[0][0] = x[bi][dim1] - box_size[dim1] / 2 + box_center[dim1];
                plane_limits[0][1] = x[bi][dim1] + box_size[dim1] / 2 + box_center[dim1];
                plane_limits[1][0] = x[bi][dim2] - box_size[dim2] / 2 + box_center[dim2];
                plane_limits[1][1] = x[bi][dim2] + box_size[dim2] / 2 + box_center[dim2];

                intersection_flag = 1;
                //compute perpendicular projection of the line connecting i and j
                proj = x[i][coordflag] - planecoord;

                //use the ratio between this projection and the i-j displacement normal to the plane
                //to define the line parameter (0-1) at the point of intersection
                lparam = proj / (delxa[coordflag]);
                if (delxa[coordflag] == 0) intersection_flag = 0;

                //use line parameter to extrapolate the possible intersection point between i-j
                intersection_point[dim1] = xj[dim1] + delxa[dim1] * (1 - lparam);
                intersection_point[dim2] = xj[dim2] + delxa[dim2] * (1 - lparam);

                //test if i-j normal coordinates are on opposing sides of the plane
                if ((proj < 0 && x[j][coordflag] < planecoord) ||
                    ((proj > 0) && x[j][coordflag] > planecoord))
                  intersection_flag = 0;

                //test the tangential coordinates to determine if the line through i-j crosses the finite sized plane
                if (intersection_point[dim1] <= plane_limits[0][0] ||
                    intersection_point[dim1] > plane_limits[0][1])
                  intersection_flag = 0;
                if (intersection_point[dim2] <= plane_limits[1][0] ||
                    intersection_point[dim2] > plane_limits[1][1])
                  intersection_flag = 0;
                //a line can only intersect a box twice; ties at edges and corners are chosen via first come first serve
                if (intersection_flag) intersect_count++;
                if (intersection_flag) {
                  if (x[i][coordflag] < planecoord) {
                    sign = -normal_flag;
                  } else {
                    sign = normal_flag;
                  }
                  atom_fluxes[bi][4 * isl] +=
                      (interaction_forceij[0] * vix + interaction_forceij[1] * viy +
                       interaction_forceij[2] * viz - interaction_forceji[0] * vjx -
                       interaction_forceji[1] * vjy - interaction_forceji[2] * vjz) *
                      sign;
                  atom_fluxes[bi][4 * isl + 1] -=
                      (interaction_forceij[0] - interaction_forceji[0]) * sign;
                  atom_fluxes[bi][4 * isl + 2] -=
                      (interaction_forceij[1] - interaction_forceji[1]) * sign;
                  atom_fluxes[bi][4 * isl + 3] -=
                      (interaction_forceij[2] - interaction_forceji[2]) * sign;
                }
                if (intersect_count == 2) break;
              }
            } else {
              //test if the bounding box (skin added) for this atom intersects with the line between i and j
              //loop over the 6 bounding planes of the box
              if (mask[bi] & groupbit) {
                intersect_count = 0;
                //determine which of the 6 planes of the atom box are intersected by this i-j pair
                for (int isl = 0; isl < 2 * domain->dimension; isl++) {
                  int is = isl / 2;
                  if (isl % 2 == 0)
                    normal_flag = 1;
                  else
                    normal_flag = -1;
                  coordflag = is;
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

                  //test negative and positive sides of the box dimension
                  if (isl % 2 == 0)
                    planecoord = x[bi][is] - box_size[is] / 2 + box_center[is];
                  else
                    planecoord = x[bi][is] + box_size[is] / 2 + box_center[is];
                  plane_limits[0][0] = x[bi][dim1] - box_size[dim1] / 2 + box_center[dim1];
                  plane_limits[0][1] = x[bi][dim1] + box_size[dim1] / 2 + box_center[dim1];
                  plane_limits[1][0] = x[bi][dim2] - box_size[dim2] / 2 + box_center[dim2];
                  plane_limits[1][1] = x[bi][dim2] + box_size[dim2] / 2 + box_center[dim2];

                  intersection_flag = 1;
                  //compute perpendicular projection of the line connecting i and j
                  proj = x[i][coordflag] - planecoord;

                  //use the ratio between this projection and the i-j displacement normal to the plane
                  //to define the line parameter (0-1) at the point of intersection
                  lparam = proj / (delxa[coordflag]);
                  if (delxa[coordflag] == 0) intersection_flag = 0;

                  //use line parameter to extrapolate the possible intersection point between i-j
                  intersection_point[dim1] = xj[dim1] + delxa[dim1] * (1 - lparam);
                  intersection_point[dim2] = xj[dim2] + delxa[dim2] * (1 - lparam);

                  //test if i-j normal coordinates are on opposing sides of the plane
                  if ((proj < 0 && x[j][coordflag] < planecoord) ||
                      ((proj > 0) && x[j][coordflag] > planecoord))
                    intersection_flag = 0;

                  //test the tangential coordinates to determine if the line through i-j crosses the finite sized plane
                  if (intersection_point[dim1] <= plane_limits[0][0] ||
                      intersection_point[dim1] > plane_limits[0][1])
                    intersection_flag = 0;
                  if (intersection_point[dim2] <= plane_limits[1][0] ||
                      intersection_point[dim2] > plane_limits[1][1])
                    intersection_flag = 0;
                  //a line can only intersect a box twice; ties at edges and corners are chosen via first come first serve
                  if (intersection_flag) intersect_count++;
                  if (intersection_flag) {
                    if (x[i][coordflag] < planecoord) {
                      sign = -normal_flag;
                    } else {
                      sign = normal_flag;
                    }
                    atom_fluxes[bi][4 * isl] +=
                        (interaction_forceij[0] * vix + interaction_forceij[1] * viy +
                         interaction_forceij[2] * viz - interaction_forceji[0] * vjx -
                         interaction_forceji[1] * vjy - interaction_forceji[2] * vjz) *
                        sign;
                    atom_fluxes[bi][4 * isl + 1] -=
                        (interaction_forceij[0] - interaction_forceji[0]) * sign;
                    atom_fluxes[bi][4 * isl + 2] -=
                        (interaction_forceij[1] - interaction_forceji[1]) * sign;
                    atom_fluxes[bi][4 * isl + 3] -=
                        (interaction_forceij[2] - interaction_forceji[2]) * sign;
                  }
                  if (intersect_count == 2) break;
                }
                //add atom box to the set of intersects for this pair if an intersection was found
                if (intersect_count) {
                  //check if max_intersect for this pair has been reached and reallocate accordingly
                  if (ij_nplane_intersects == max_intersect_count) {
                    max_intersect_count += MAXBINOVERLAP;
                    memory->grow(ij_plane_intersects, max_intersect_count,
                                 "FixSurfaceFluxesAtom:ij_plane_intersects");
                  }
                  ij_plane_intersects[ij_nplane_intersects++] = bi;
                  scan_flags[bi] = 1;
                }
              }
            }
          }

          //consider if more bin intersections are possible
          if (grid_travel[0] == -delta_grid[0] && grid_travel[1] == -delta_grid[1] &&
              grid_travel[2] == -delta_grid[2])
            break;

          //test bin surface planes for intersection with the connecting line and move to another bin accordingly
          intersection_flags[0] = intersection_flags[1] = intersection_flags[2] = 0;

          //an intersection must take place
          for (int is = 0; is < domain->dimension; is++) {
            if (grid_travel[is] == -delta_grid[is]) continue;
            coordflag = is;

            //displace grid indices by 1 if testing the negative direction
            if (delta_grid[is] > 0) ix[is]--;

            grid_point[0] = binsize[0] * ix[0] + bboxlo[0];
            grid_point[1] = binsize[1] * ix[1] + bboxlo[1];
            grid_point[2] = binsize[2] * ix[2] + bboxlo[2];
            planecoord = grid_point[is] + binsize[is];
            plane_limits[0][0] = grid_point[dim1];
            plane_limits[0][1] = grid_point[dim1] + binsize[dim1];
            plane_limits[1][0] = grid_point[dim2];
            plane_limits[1][1] = grid_point[dim2] + binsize[dim2];

            //intersection_flags[is] = 1;
            //compute perpendicular projection of the line connecting i and j
            proj = x[i][coordflag] - planecoord;

            //use the ratio between this projection and the i-j displacement normal to the plane
            //to define the line parameter (0-1) at the point of intersection\

            if (delxa[is] < DELEPSILON && delxa[is] > -DELEPSILON) dels[is] = 0;
            dels[is] = proj / (delxa[coordflag]);

            //reset ix
            if (delta_grid[is] > 0) ix[is]++;
          }

          //set intersection flag to be the direction with the minimum line traversal
          //initialize to the first direction that still has to travel
          if (grid_travel[0] != -delta_grid[0]) {
            smin = dels[0];
            mindir = 0;
          } else if (grid_travel[1] != -delta_grid[1]) {
            smin = dels[1];
            mindir = 1;
          } else {
            smin = dels[2];
            mindir = 2;
          }

          //compare to remaining directions that still need to travel and decide
          for (int is = mindir + 1; is < domain->dimension; is++) {
            if (grid_travel[is] == -delta_grid[is]) continue;

            if (dels[is] < smin) {
              smin = dels[is];
              mindir = is;
            }
          }
          intersection_flags[mindir] = 1;

          //if intersection flag was set. Only one direction at a time is allowed.
          if (intersection_flags[0]) {
            if (delta_grid[0] < 0)
              grid_travel[0]++;
            else
              grid_travel[0]--;
          } else if (intersection_flags[1]) {
            if (delta_grid[1] < 0)
              grid_travel[1]++;
            else
              grid_travel[1]--;
          } else if (intersection_flags[2]) {
            if (delta_grid[2] < 0)
              grid_travel[2]++;
            else
              grid_travel[2]--;
          }

          if (iloop == max_intersection)
            error->one(FLERR, "grid search did not reach end point; fix the code.");
        }

        //reset scan flags for the next i-j pair
        for (int reset = 0; reset < ij_nplane_intersects; reset++)
          scan_flags[ij_plane_intersects[reset]] = 0;
      }
    }
  }

  comm->reverse_comm_fix(this);

  //convert units of force to stress times area units
  double nktv2p = force->nktv2p;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      for (j = 0; j < 24; j++) {
        if (j % 4 != 0) atom_fluxes[i][j] *= nktv2p;
      }
    }
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesAtomOp::box2bins(int atom_index)
{
  double bounding_boxlo[3];
  double bounding_boxhi[3];
  double skin = neighbor->skin;
  double *cutghost = comm->cutghost;
  double *x = atom->x[atom_index];
  int ixl, iyl, izl, ixh, iyh, izh;

  //calculate the bounding box for this atom

  bounding_boxlo[0] = x[0] - box_size[0] / 2 + box_center[0];
  bounding_boxlo[1] = x[1] - box_size[1] / 2 + box_center[1];
  bounding_boxlo[2] = x[2] - box_size[2] / 2 + box_center[2];
  bounding_boxhi[0] = x[0] + box_size[0] / 2 + box_center[0];
  bounding_boxhi[1] = x[1] + box_size[1] / 2 + box_center[1];
  bounding_boxhi[2] = x[2] + box_size[2] / 2 + box_center[2];

  //compute lowest bin overlap id for each dimension
  if (bounding_boxlo[0] >= bboxhi[0])
    ixl = static_cast<int>((bounding_boxlo[0] - bboxhi[0]) * bininvx) + nbinx;
  else if (bounding_boxlo[0] >= bboxlo[0]) {
    ixl = static_cast<int>((bounding_boxlo[0] - bboxlo[0]) * bininvx);
    ixl = min(ixl, nbinx - 1);
  } else
    ixl = static_cast<int>((bounding_boxlo[0] - bboxlo[0]) * bininvx) - 1;

  if (bounding_boxlo[1] >= bboxhi[1])
    iyl = static_cast<int>((bounding_boxlo[1] - bboxhi[1]) * bininvy) + nbiny;
  else if (bounding_boxlo[1] >= bboxlo[1]) {
    iyl = static_cast<int>((bounding_boxlo[1] - bboxlo[1]) * bininvy);
    iyl = min(iyl, nbiny - 1);
  } else
    iyl = static_cast<int>((bounding_boxlo[1] - bboxlo[1]) * bininvy) - 1;

  if (bounding_boxlo[2] >= bboxhi[2])
    izl = static_cast<int>((bounding_boxlo[2] - bboxhi[2]) * bininvz) + nbinz;
  else if (bounding_boxlo[2] >= bboxlo[2]) {
    izl = static_cast<int>((bounding_boxlo[2] - bboxlo[2]) * bininvz);
    izl = min(izl, nbinz - 1);
  } else
    izl = static_cast<int>((bounding_boxlo[2] - bboxlo[2]) * bininvz) - 1;

  //compute highest bin overlap id for each dimension
  if (bounding_boxhi[0] >= bboxhi[0])
    ixh = static_cast<int>((bounding_boxhi[0] - bboxhi[0]) * bininvx) + nbinx;
  else if (bounding_boxhi[0] >= bboxlo[0]) {
    ixh = static_cast<int>((bounding_boxhi[0] - bboxlo[0]) * bininvx);
    ixh = min(ixh, nbinx - 1);
  } else
    ixh = static_cast<int>((bounding_boxhi[0] - bboxlo[0]) * bininvx) - 1;

  if (bounding_boxhi[1] >= bboxhi[1])
    iyh = static_cast<int>((bounding_boxhi[1] - bboxhi[1]) * bininvy) + nbiny;
  else if (bounding_boxhi[1] >= bboxlo[1]) {
    iyh = static_cast<int>((bounding_boxhi[1] - bboxlo[1]) * bininvy);
    iyh = min(iyh, nbiny - 1);
  } else
    iyh = static_cast<int>((bounding_boxhi[1] - bboxlo[1]) * bininvy) - 1;

  if (bounding_boxhi[2] >= bboxhi[2])
    izh = static_cast<int>((bounding_boxhi[2] - bboxhi[2]) * bininvz) + nbinz;
  else if (bounding_boxhi[2] >= bboxlo[2]) {
    izh = static_cast<int>((bounding_boxhi[2] - bboxlo[2]) * bininvz);
    izh = min(izh, nbinz - 1);
  } else
    izh = static_cast<int>((bounding_boxhi[2] - bboxlo[2]) * bininvz) - 1;

  //binning bounding boxes of ghosts can produce bin indices exceeding the
  //local bin limits; set overlap to the limit in this case. Note that these
  //limiting cases do not contribute to the flux calculation anyway since ghosts
  //to be looped over lie within the bin limits.
  if (ixl < mbinxlo) ixl = mbinxlo;
  if (iyl < mbinylo) iyl = mbinylo;
  if (izl < mbinzlo) izl = mbinzlo;
  if (ixh > mbinxhi) ixh = mbinxhi;
  if (iyh > mbinyhi) iyh = mbinyhi;
  if (izh > mbinzhi) izh = mbinzhi;

  bin_overlap_limits[0] = ixl;
  bin_overlap_limits[1] = iyl;
  bin_overlap_limits[2] = izl;
  bin_overlap_limits[3] = ixh;
  bin_overlap_limits[4] = iyh;
  bin_overlap_limits[5] = izh;
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesAtomOp::allocate_arrays(int n, int mode)
{
  //grows bin arrays
  if (mode == 0) {
    if (!first_alloc) {
      first_alloc = 1;
      memory->grow(bin_expansion_counts, n, "FixSurfaceFluxesAtom:bin_ncontent");
      for (int init = 0; init < n; init++) { bin_expansion_counts[init] = max_bin_expansion_count; }
      bin_content = (int **) memory->smalloc(n * sizeof(int *), "FixSurfaceFluxesAtom:bin_content");
      for (int init = 0; init < n; init++)
        memory->create(bin_content[init], MAXBINCONTENT, "FixSurfaceFluxesAtom:bin_content");
    } else {
      //compute maximum expansion count used in previous allocations
      for (int init = 0; init < maxbin; init++) {
        if (bin_expansion_counts[init] > max_bin_expansion_count)
          max_bin_expansion_count = bin_expansion_counts[init];
      }
      bin_content = (int **) memory->srealloc(bin_content, n * sizeof(int *),
                                              "FixSurfaceFluxesAtom:bin_content");
      for (int init = 0; init < n; init++) {
        if (init >= maxbin) bin_content[init] = NULL;
        memory->grow(bin_content[init], MAXBINCONTENT + max_bin_expansion_count * EXPAND,
                     "FixSurfaceFluxesAtom:bin_content");
      }
      memory->grow(bin_expansion_counts, n, "FixSurfaceFluxesAtom:bin_ncontent");
      for (int init = 0; init < n; init++) bin_expansion_counts[init] = max_bin_expansion_count;
    }
    maxbin = n;
    memory->grow(bin_ncontent, maxbin, "FixSurfaceFluxesAtom:bin_ncontent");
  }
}

/* ---------------------------------------------------------------------- */

int FixSurfaceFluxesAtomOp::pack_reverse_comm(int n, int first, double *buf)
{
  int i, l, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    for (l = 0; l < 24; l++) buf[m++] = atom_fluxes[i][l];
  return m;
}

/* ---------------------------------------------------------------------- */

void FixSurfaceFluxesAtomOp::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, l, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (l = 0; l < 24; l++) atom_fluxes[j][l] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double FixSurfaceFluxesAtomOp::memory_usage()
{
  double bytes = 0;

  bytes += memory->usage(atom_fluxes, atom->nlocal + atom->nghost, 24);
  return bytes;
}