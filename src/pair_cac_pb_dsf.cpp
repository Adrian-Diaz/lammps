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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "pair_cac_pb_dsf.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "update.h"
#include "neigh_list.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "asa_user.h"

#define MAXNEIGHIN  100
#define MAXNEIGHOUT  200
#define EXPAND 10
#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429
//#include "math_extra.h"
using namespace LAMMPS_NS;
using namespace MathConst;


/* ---------------------------------------------------------------------- */

PairCACPbDSF::PairCACPbDSF(LAMMPS *lmp) : PairCAC(lmp)
{
  restartinfo = 0;
  nmax = 0;
  outer_neighflag = 0;
  flux_enable = 1;
}

/* ---------------------------------------------------------------------- */

PairCACPbDSF::~PairCACPbDSF() {
  if (allocated) {
  memory->destroy(setflag);
  memory->destroy(cutsq);
	memory->destroy(cut);
	memory->destroy(a);
	memory->destroy(rho);
	memory->destroy(c);
	memory->destroy(rhoinv);
	memory->destroy(buck1);
	memory->destroy(buck2);
	memory->destroy(offset);
	memory->destroy(inner_neighbor_coords);
	memory->destroy(inner_neighbor_types);
	memory->destroy(inner_neighbor_charges);
  }
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairCACPbDSF::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  max_nodes_per_element = atom->nodes_per_element;
  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut, n + 1, n + 1, "pair:cut_lj");
  memory->create(a, n + 1, n + 1, "pair:a");
  memory->create(rho, n + 1, n + 1, "pair:rho");
  memory->create(c, n + 1, n + 1, "pair:c");
  memory->create(rhoinv, n + 1, n + 1, "pair:rhoinv");
  memory->create(buck1, n + 1, n + 1, "pair:buck1");
  memory->create(buck2, n + 1, n + 1, "pair:buck2");
  memory->create(offset, n + 1, n + 1, "pair:offset");


  memory->create(mass_matrix,max_nodes_per_element, max_nodes_per_element,"pairCAC:mass_matrix");
  memory->create(mass_copy, max_nodes_per_element, max_nodes_per_element,"pairCAC:copy_mass_matrix");
  memory->create(force_column, max_nodes_per_element,3,"pairCAC:force_residue");
  memory->create(current_force_column, max_nodes_per_element,"pairCAC:current_force_residue");
  memory->create(current_nodal_forces, max_nodes_per_element,"pairCAC:current_nodal_force");
  memory->create(pivot, max_nodes_per_element+1,"pairCAC:pivots");
  quadrature_init(2);
}

/* ----------------------------------------------------------------------
global settings
------------------------------------------------------------------------- */
void PairCACPbDSF::settings(int narg, char **arg) {
	if (narg <3 || narg>4) error->all(FLERR, "Illegal pair_style command");

	force->newton_pair = 0;
	cut_buck = utils::numeric(FLERR, arg[0],false,lmp);
    alf = utils::numeric(FLERR, arg[1],false,lmp);
    cut_coul = utils::numeric(FLERR, arg[2],false,lmp);
	 if (narg == 4) {
		if (strcmp(arg[3], "one") == 0) atom->one_layer_flag=one_layer_flag = 1;
		else error->all(FLERR, "Unexpected argument in CAC/Pb pair style invocation");

	}
	cut_global_s = cut_buck;
	if(cut_coul>cut_buck)
	 cut_global_s=cut_coul;

	if (allocated) {
		int i, j;
		for (i = 1; i <= atom->ntypes; i++)
			for (j = i; j <= atom->ntypes; j++)
				if (setflag[i][j]) cut[i][j] = cut_global_s;
	}
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairCACPbDSF::coeff(int narg, char **arg) {
	if (narg < 5 || narg > 6)
		error->all(FLERR, "Incorrect args for pair coefficients");
	if (!allocated) allocate();
	setflag[1][1] = 1;
	int ilo, ihi, jlo, jhi;
	utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
    utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

	double a_one = utils::numeric(FLERR, arg[2],false,lmp);
  double rho_one = utils::numeric(FLERR, arg[3],false,lmp);
  if (rho_one <= 0) error->all(FLERR, "Incorrect args for pair coefficients");
  double c_one = utils::numeric(FLERR, arg[4],false,lmp);

	double cut_one = cut_global_s;
	if (narg == 6) cut_one = utils::numeric(FLERR, arg[5],false,lmp);

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo, i); j <= jhi; j++) {
			a[i][j] = a_one;
			rho[i][j] = rho_one;
			c[i][j] = c_one;
			cut[i][j] = cut_one;
			setflag[i][j] = 1;
			count++;
		}
	}

	if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairCACPbDSF::init_one(int i, int j) {

	//if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

	rhoinv[i][j] = 1.0 / rho[i][j];
	buck1[i][j] = a[i][j] / rho[i][j];
	buck2[i][j] = 6.0*c[i][j];

	if (offset_flag && (cut[i][j] > 0.0)) {
		double rexp = exp(-cut[i][j] / rho[i][j]);
		offset[i][j] = a[i][j] * rexp - c[i][j] / pow(cut[i][j], 6.0);
	}
	else offset[i][j] = 0.0;

	a[j][i] = a[i][j];
	c[j][i] = c[i][j];
	rhoinv[j][i] = rhoinv[i][j];
	buck1[j][i] = buck1[i][j];
	buck2[j][i] = buck2[i][j];
	offset[j][i] = offset[i][j];

	// compute I,J contribution to long-range tail correction
	// count total # of atoms of type I and J via Allreduce

	if (tail_flag) {
		int *type = atom->type;
		int nlocal = atom->nlocal;

		double count[2], all[2];
		count[0] = count[1] = 0.0;
		for (int k = 0; k < nlocal; k++) {
			if (type[k] == i) count[0] += 1.0;
			if (type[k] == j) count[1] += 1.0;
		}
		MPI_Allreduce(count, all, 2, MPI_DOUBLE, MPI_SUM, world);

		double rho1 = rho[i][j];
		double rho2 = rho1*rho1;
		double rho3 = rho2*rho1;
		double rc = cut[i][j];
		double rc2 = rc*rc;
		double rc3 = rc2*rc;
		etail_ij = 2.0*MY_PI*all[0] * all[1] *
			(a[i][j] * exp(-rc / rho1)*rho1*(rc2 + 2.0*rho1*rc + 2.0*rho2) -
				c[i][j] / (3.0*rc3));
		ptail_ij = (-1 / 3.0)*2.0*MY_PI*all[0] * all[1] *
			(-a[i][j] * exp(-rc / rho1)*
			(rc3 + 3.0*rho1*rc2 + 6.0*rho2*rc + 6.0*rho3) + 2.0*c[i][j] / rc3);
	}


	return cut_global_s;
}

/* ---------------------------------------------------------------------- */


void PairCACPbDSF::init_style()
{
	PairCAC::init_style();
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style CAC_Buck requires atom IDs");
  if (!atom->q_flag)
	  error->all(FLERR, "Pair coul/wolf requires atom attribute q for charges");
  atom->max_neigh_inner_init = maxneigh_quad_inner = MAXNEIGHIN;
  atom->max_neigh_outer_init = maxneigh_quad_outer = MAXNEIGHOUT;
  cut_coulsq = cut_coul*cut_coul;
  double erfcc = erfc(alf*cut_coul);
  double erfcd = exp(-alf*alf*cut_coul*cut_coul);
  f_shift = -(erfcc/cut_coulsq + 2.0/MY_PIS*alf*erfcd/cut_coul);
  e_shift = erfcc/cut_coul - f_shift*cut_coul;
}

//-----------------------------------------------------------------------

void PairCACPbDSF::force_densities(int iii, double s, double t, double w, double coefficients,
	double &force_densityx, double &force_densityy, double &force_densityz) {

int internal;
double delx,dely,delz;
double r2inv;
double r6inv;
double ecoul;
int timestep=update->ntimestep;
double *special_lj = force->special_lj;
double  forcebuck, factor_lj, fpair;
double  fpair2, forcecoul, factor_coul;
double prefactor;
double r, rexp;
int *type = atom->type;
double distancesq;
double scan_position[3];
double rcut;
int current_type = poly_counter;
double erfcc, erfcd, v_sh, dvdrr, e_self, e_shift, f_shift, qisq;
double *special_coul = force->special_coul;
double qqrd2e = force->qqrd2e;

int flagm;
int nodes_per_element;

	rcut = cut_global_s;
	int origin_type = type_array[poly_counter];
	int listtype;
	int listindex;
	int poly_index;
	int scan_type;
	int element_index;
	int neigh_max = inner_quad_lists_counts[pqi];
	double ****nodal_positions = atom->nodal_positions;
	int **node_types = atom->node_types;
	double **node_charges = atom->node_charges;
	double origin_element_charge = node_charges[iii][poly_counter];
	double neighbor_element_charge;
  int **inner_quad_indices = inner_quad_lists_index[pqi];
	qisq = origin_element_charge*origin_element_charge;

	  //allocate arrays that store neighbor information around just this quadrature point
      allocate_quad_memory();
      //set virtual neighbor types, etc.
      init_quad_arrays();
      //interpolate virtual atom coordinates from shape functions corresponding to unit cells
      interpolation(iii,s,t,w);
			for (int l = 0; l < neigh_max; l++) {

				scan_type = inner_neighbor_types[l];
				scan_position[0] = inner_neighbor_coords[l][0];
				scan_position[1] = inner_neighbor_coords[l][1];
				scan_position[2] = inner_neighbor_coords[l][2];
				neighbor_element_charge = inner_neighbor_charges[l];

				delx = current_position[0] - scan_position[0];
				dely = current_position[1] - scan_position[1];
				delz = current_position[2] - scan_position[2];
				distancesq = delx*delx + dely*dely + delz*delz;
				if (distancesq < cut_global_s*cut_global_s) {
					fpair = 0;
					fpair2 = 0;
					if (distancesq < cut_coul*cut_coul) {
						factor_coul = special_coul[sbmask(inner_quad_indices[l][0])];
						r = sqrt(distancesq);
						prefactor = qqrd2e*origin_element_charge*neighbor_element_charge / r;
                        erfcd = exp(-alf*alf*distancesq);
                        t = 1.0 / (1.0 + EWALD_P*alf*r);
                        erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd;
                        forcecoul = prefactor * (erfcc/r + 2.0*alf/MY_PIS * erfcd +
                          r*f_shift) * r;
                        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
                          fpair = forcecoul / distancesq;
                        if (quad_eflag){
                          ecoul = prefactor * (erfcc - r*e_shift - distancesq*f_shift);
                        if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
            }
					}
					if (distancesq < cut_buck*cut_buck) {
						if ((scan_type==1&&origin_type!=1)|| (scan_type != 1 && origin_type == 1)|| (scan_type != 1 && origin_type != 1)) {
							//buckingham force
							r2inv = 1.0 / distancesq;
							r6inv = r2inv*r2inv*r2inv;
							rexp = exp(-r*rhoinv[origin_type][scan_type]);
							forcebuck = buck1[origin_type][scan_type] * r*rexp - buck2[origin_type][scan_type] * r6inv;
							fpair2 = forcebuck*r2inv;
						}
					}
					//compute total pair force
					force_densityx += delx*(fpair + fpair2);
					force_densityy += dely*(fpair + fpair2);
					force_densityz += delz*(fpair + fpair2);
					if(atom->CAC_virial){
		      virial_density[0] += 0.5*delx*delx*fpair;
		      virial_density[1] += 0.5*dely*dely*fpair;
		      virial_density[2] += 0.5*delz*delz*fpair;
		      virial_density[3] += 0.5*delx*dely*fpair;
		      virial_density[4] += 0.5*delx*delz*fpair;
		      virial_density[5] += 0.5*dely*delz*fpair;
		      }
          if (quad_eflag)
          quadrature_energy += (a[origin_type][scan_type]*rexp - c[origin_type][scan_type]*r6inv -
            offset[origin_type][scan_type])/2+ecoul/2;
				}
			}

  //additional cac flux contributions due to neighbors interacting with neighbors
  //  in the vicinity of this quadrature point
  if (quad_flux_flag) {
    //compute_intersections();
    quad_neigh_flux();
  }

}

/* ---------------------------------------------------------------------- */

double PairCACPbDSF::pair_interaction_q(double distancesq, int itype, int jtype
                                          , double qi, double qj)
{
  double fpair, prefactor, dvdrr;
  double r, erfcc, erfcd, t;
  double qqrd2e = force->qqrd2e;
  double forcebuck = 0;
  double forcecoul = 0;
  double r2inv;
  double v_sh, rexp, r6inv;

  r = sqrt(distancesq);
  if(distancesq < cut_buck*cut_buck){
  if ((jtype==1&&itype!=1)|| (jtype != 1 && itype == 1)|| (jtype != 1 && itype != 1)) {
  r2inv = 1.0 / distancesq;
  r6inv = r2inv*r2inv*r2inv;
  rexp = exp(-r*rhoinv[itype][jtype]);
  forcebuck = buck1[itype][jtype] * r*rexp - buck2[itype][jtype] * r6inv;
  forcebuck = forcebuck*r2inv;
  }
  }
  if(distancesq < cut_coul*cut_coul){
  prefactor = qqrd2e*qi*qj / r;
  erfcd = exp(-alf*alf*distancesq);
  t = 1.0 / (1.0 + EWALD_P*alf*r);
  erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd;
  forcecoul = prefactor * (erfcc/r + 2.0*alf/MY_PIS * erfcd + r*f_shift) * r;
  //if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul)*prefactor;
  forcecoul = forcecoul / distancesq;
  //if(quad_eflag)
    //if (factor_coul < 1.0) v_sh -= (1.0-factor_coul)*prefactor;
  }
  fpair = forcecoul + forcebuck;
  return fpair;
}
