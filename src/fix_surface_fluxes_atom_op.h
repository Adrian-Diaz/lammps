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

FixStyle(surface/fluxes/atom/op,FixSurfaceFluxesAtomOp)

#else

#ifndef LMP_FIX_SURFACE_FLUXES_ATOM_OP_H
#define LMP_FIX_SURFACE_FLUXES_ATOM_OP_H

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixSurfaceFluxesAtomOp : public Fix {
 public:
  FixSurfaceFluxesAtomOp(class LAMMPS *, int, char **);
  ~FixSurfaceFluxesAtomOp();

  
  double **atom_fluxes; //the fluxes corresponding to each atoms surrounding box
  int *bin_ncontent;          //stores how many atom boxes overlapped with each bin
  int **bin_content;          //stores which atom boxes overlap with each bin
  int *ij_plane_intersects; //stores which atom box planes the current bond line intersects
  int ij_nplane_intersects; //stores how many atom box planes the current bond line intersects

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
  double binsizex,binsizey,binsizez,binsize[3];  // bin sizes and inverse sizes
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

*/
