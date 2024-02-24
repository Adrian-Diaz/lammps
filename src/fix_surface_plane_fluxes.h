/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(surface/plane/fluxes,FixSurfacePlaneFluxes)

#else

#ifndef LMP_FIX_SURFACE_PLANE_FLUXES_H
#define LMP_FIX_SURFACE_PLANE_FLUXES_H

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixSurfacePlaneFluxes : public Fix {
 public:
  FixSurfacePlaneFluxes(class LAMMPS *, int, char **);
  ~FixSurfacePlaneFluxes();
  int setmask();
  void init();

  virtual void post_neighbor();
  virtual void setup_post_neighbor();
  virtual double compute_vector(int);


 private:
  //bin data
  class NBin *nb;               // ptr to bin object
  int nbinx,nbiny,nbinz;           // # of global bins
  int mbins;                       // # of local bins and offset on this proc
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;
  int mbinxhi,mbinyhi,mbinzhi;
  double binsizex,binsizey,binsizez;  // bin sizes and inverse sizes
  double bininvx,bininvy,bininvz;
  int *binhead;                // index of first atom in each bin
  int *bins;                   // index of next atom in same bin
  int *atom2bin;               // bin assignment for each atom (local+ghost)
  int bin_limits[3][2];     //bin extents to loop over for a single plane
  int ****grid_bin_limits;   //bin extents to loop over for every plane in the grid
  int invoked;

  //stencil data
  class NStencil *ns;           // ptr to NStencil instance I depend on
  int nstencil;
  int *stencil;
  int **stencilxyz;
  int *nstencil_multi;
  int **stencil_multi;
  double **distsq_multi;
  int sx,sy,sz;

  //printing related
  FILE *fp;
  int pinterval;
  int tecplot_flag;   //denotes whether the output format is Tecplot points
  int ntimestep,ptimestep;
  int surface_select;  //denotes which grid surface normals are printed if user specifies preference
  
  int local_computation;
  int irequest;         //neighbor request index for this fix
  int halforfull;       //0 if neighbor list being used is half and 1 if full
  int avg_interval, avg_divisor;
  int dimension;
  char *group2;
  char* cordflag;
  char *id_pe;
  int me;
  int box_count[3];
  double grid_size[3];
  double *bboxlo;
  double *bboxhi;
  int coordflag;
  int gridflag;
  int sub_grid_limits[3][2];
  double tangentflag;
  double surface_size_x,surface_size_y;
  int surface_count_x;
  int surface_count_y;
  class Pair *pair;
  double ***surface_sums_2D; //2D grid of surface elements
  double **surface_sums_1D;  //1D grid of surface elements
  double ***surface_sums_2D_global; //2D grid of surface elements
  double **surface_sums_1D_global;  //1D grid of surface elements
  int planecount;
  int jgroup,jgroupbit,othergroupbit;
  int nmax;
  class Compute *ccpe;
  double **cutsq;
  double planecoord;
  double *interaction_force;
  double *eatom;
  double e_self,e_correction;
  int pairflag,kspaceflag,boundaryflag;
  int columncount;
  double *fluxes;
  double *cfluxes;
  double ***grid_surface_fluxes;
  double ***global_grid_surface_fluxes;
  int **grid_compute_flags;
  int *gather_counts;
  double debug_sum;
  double plane_limits[2][2];

  void bin_info();
  void pair_contribution();
  inline int sbmask(int j) const {
    return j >> SBBITS & 3;
  }

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: No values in fix ave/time command

Self-explanatory.

E: Compute ID for fix ave/time does not exist

Self-explanatory.

E: Fix ID for fix ave/time does not exist

Self-explanatory.

E: Invalid fix ave/time off column

Self-explantory.

E: Fix ave/time compute does not calculate a scalar

Self-explantory.

E: Fix ave/time compute does not calculate a vector

Self-explantory.

E: Fix ave/time compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix ave/time compute does not calculate an array

Self-explanatory.

E: Fix ave/time compute array is accessed out-of-range

An index for the array is out of bounds.

E: Fix ave/time fix does not calculate a scalar

Self-explanatory.

E: Fix ave/time fix does not calculate a vector

Self-explanatory.

E: Fix ave/time fix vector cannot be variable length

Self-explanatory.

E: Fix ave/time fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix for fix ave/time not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/time
is requesting a value on a non-allowed timestep.

E: Fix ave/time fix does not calculate an array

Self-explanatory.

E: Fix ave/time fix array cannot be variable length

Self-explanatory.

E: Fix ave/time fix array is accessed out-of-range

An index for the array is out of bounds.

E: Variable name for fix ave/time does not exist

Self-explanatory.

E: Fix ave/time variable is not equal-style variable

Self-explanatory.

E: Fix ave/time cannot use variable with vector mode

Variables produce scalar values.

E: Error writing file header

Something in the output to the file triggered an error.

E: Fix ave/time cannot set output array intensive/extensive from these inputs

One of more of the vector inputs has individual elements which are
flagged as intensive or extensive.  Such an input cannot be flagged as
all intensive/extensive when turned into an array by fix ave/time.

E: Invalid timestep reset for fix ave/time

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

E: Error writing out time averaged data

Something in the output to the file triggered an error.

E: Fix ave/time columns are inconsistent lengths

Self-explanatory.

E: Cannot open fix ave/time file %s

The specified file cannot be opened.  Check that the path and name are
correct.

*/
