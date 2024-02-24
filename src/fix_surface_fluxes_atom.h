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

FixStyle(surface/fluxes/atom,FixSurfaceFluxesAtom)

#else

#ifndef LMP_FIX_SURFACE_FLUXES_ATOM_H
#define LMP_FIX_SURFACE_FLUXES_ATOM_H

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixSurfaceFluxesAtom : public Fix {
 public:
  FixSurfaceFluxesAtom(class LAMMPS *, int, char **);
  ~FixSurfaceFluxesAtom();

  
  double **atom_fluxes; //the fluxes corresponding to each atoms surrounding box
  int *bin_ncontent;          //stores how many atom boxes overlapped with each bin
  int **bin_content;          //stores which atom boxes overlap with each bin
  int ***ij_plane_intersects; //stores which atom box planes the current bond line intersects
  int **ij_nplane_intersects; //stores how many atom box planes the current bond line intersects

  int setmask();
  void init();
  virtual void end_of_step();
  virtual void post_neighbor();
  virtual void setup_post_neighbor();
  virtual void setup(int);
  virtual int pack_reverse_comm(int, int, double *);
  virtual void unpack_reverse_comm(int, int *, double *);
  virtual double memory_usage();

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
  double *bboxlo,*bboxhi;
  int *binhead;                // index of first atom in each bin
  int *bins;                   // index of next atom in same bin
  int *atom2bin;               // bin assignment for each atom (local+ghost)
  int *group_include;
  int *bin_overlap_limits;
  int first_alloc;
  int maxbin;
  int max_bin_expansion_count;
  int *bin_expansion_counts;
  int *max_neighbor_counts;
  int *old_max_neighbor_counts;
  int **max_intersect_counts;
  int max_intersect_count;

  //stencil data
  class NStencil *ns;           // ptr to NStencil instance I depend on
  int nstencil;
  int *stencil;
  int **stencilxyz;
  int *nstencil_multi;
  int **stencil_multi;
  double **distsq_multi;
  int sx,sy,sz;

  int irequest;         //neighbor request index used by this fix
  int avg_interval, avg_divisor;
  int dimension;
  char *group2;
  char* cordflag;
  char *id_pe;
  int me;
  double box_size[3], box_center[3];
  double tangentflag;
  class Pair *pair;
  int planecount;
  int jgroup,jgroupbit,othergroupbit;
  int nmax;
  class Compute *ccpe;
  double **cutsq;
  double *interaction_force;
  double *eatom;
  double e_self,e_correction;
  int pairflag,kspaceflag,boundaryflag;
  int columncount;
  double **fluxes;
  double *cfluxes;
  int *gather_counts;
  int *scan_flags;         //used to find repeats in the intersection list due to bin search redundancies
  double debug_sum;

  void bin_info();
  void group_inclusion();
  void compute_flux();
  void box2bins(int);
  void allocate_arrays(int,int);
  int min(int a, int b){ 
    if(a<b) return a;
    else return b;
    }
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
