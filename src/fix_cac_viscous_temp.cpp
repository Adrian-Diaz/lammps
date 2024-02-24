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
The command controls the temperature of the atoms to ensure that it does 
not exceed a designated value.  (created by Yang Li)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "fix_cac_viscous_temp.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"

#include "group.h"
#include "compute.h"
#include "modify.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixViscousCACTemp::FixViscousCACTemp(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5) error->all(FLERR,"Illegal fix cac/viscous/temp command");

  double gamma_one = utils::numeric(FLERR,arg[3],false,lmp);
  gamma = new double[atom->ntypes+1];
  for (int i = 1; i <= atom->ntypes; i++) gamma[i] = gamma_one;

  t_stop = utils::numeric(FLERR,arg[4],false,lmp);
  t_current= 0;
  
  respa_level_support = 1;
  ilevel_respa = 0;

  nevery = utils::inumeric(FLERR,arg[5],false,lmp);
  if (nevery <= 0) error->all(FLERR,"Illegal fix cac/visc/temp command");
  
  // create a new compute temp
  // id = fix-ID + temp, compute group = fix group

  int n = strlen(id) + 16;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_cac_nodal_temp");

  char **newarg = new char*[6];
  newarg[0] = id_temp;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "cac/nodal/temp";
  modify->add_compute(3,newarg);
  delete [] newarg;
  tflag = 1;
}

/* ---------------------------------------------------------------------- */

FixViscousCACTemp::~FixViscousCACTemp()
{
  delete [] gamma;
  // delete temperature if fix created it

  if (tflag) modify->delete_compute(id_temp);
  delete [] id_temp;
}

/* ---------------------------------------------------------------------- */

int FixViscousCACTemp::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixViscousCACTemp::init()
{
  int max_respa = 0;
  if (!atom->CAC_flag) error->all(FLERR,"CAC fix styles require a CAC atom style");
  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = max_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,max_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixViscousCACTemp::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixViscousCACTemp::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixViscousCACTemp::post_force(int vflag)
{
  if (update->ntimestep > (update->ntimestep/nevery)*nevery)
	t_current = 0;
  else {
  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix cac/viscous/temp does not exist");
  temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,
                 "Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Group for fix_modify temp != fix group");
  
    t_current = temperature->compute_scalar();
  }
   
  // apply drag force to atoms in group
  // direction is opposed to velocity vector
  // magnitude depends on atom type

  double ****nodal_velocities = atom->nodal_velocities;
  double ****nodal_forces = atom->nodal_forces;
  int *nodes_count_list = atom->nodes_per_element_list;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *poly_count = atom->poly_count;
  int *element_type = atom->element_type;
  int **node_types = atom->node_types;
  int **element_scale = atom->element_scale;
  int nodes_per_element;
  double drag;
  
  if ( t_current > t_stop ) {
  
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      nodes_per_element = nodes_count_list[element_type[i]];

    for (int k = 0; k < poly_count[i]; k++) {
      drag = gamma[node_types[i][k]];

      for (int j = 0; j < nodes_per_element; j++) {
        nodal_forces[i][k][j][0] -= drag * nodal_velocities[i][k][j][0];
        nodal_forces[i][k][j][1] -= drag * nodal_velocities[i][k][j][1];
        nodal_forces[i][k][j][2] -= drag * nodal_velocities[i][k][j][2];
      }
    }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixViscousCACTemp::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixViscousCACTemp::min_post_force(int vflag)
{
  post_force(vflag);
}
