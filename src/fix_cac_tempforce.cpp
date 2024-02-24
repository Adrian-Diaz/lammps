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
  The command implements the thermal expansion of CAC CG elements (created by Yang Li)
------------------------------------------------------------------------- */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include "fix_cac_tempforce.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixCAC_Temp_Force::FixCAC_Temp_Force(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xstr(NULL), ystr(NULL), zstr(NULL), idregion(NULL), sforce(NULL)
{
  if (narg < 9) error->all(FLERR,"Illegal fix cac/tempforce command");

  dynamic_group_allow = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  xstr = ystr = zstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else if (strcmp(arg[3],"NULL") == 0) {
    xstyle = NONE;
  } else {
    xvalue = utils::numeric(FLERR,arg[3],false,lmp);
    xstyle = CONSTANT;
  }
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else if (strcmp(arg[4],"NULL") == 0) {
    ystyle = NONE;
  } else {
    yvalue = utils::numeric(FLERR,arg[4],false,lmp);
    ystyle = CONSTANT;
  }
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else if (strcmp(arg[5],"NULL") == 0) {
    zstyle = NONE;
  } else {
    zvalue = utils::numeric(FLERR,arg[5],false,lmp);
    zstyle = CONSTANT;
  }
  
  //element_size = utils::inumeric(FLERR,arg[6],false,lmp);
  //if (element_size <= 0) error->all(FLERR,"Illegal element_size");

  unitcell_volume = utils::numeric(FLERR,arg[6],false,lmp);
  if (unitcell_volume <= 0) error->all(FLERR,"Illegal unitcell_volume");
  
  target_temp = utils::numeric(FLERR,arg[7],false,lmp);
  if (target_temp <= 0) error->all(FLERR,"Illegal target tempurature");
 
  amp_factor = utils::numeric(FLERR,arg[8],false,lmp);
  if (amp_factor <= 0) error->all(FLERR,"Illegal amp_factor");
  // optional args

  iregion = -1;
  idregion = NULL;

  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix cac/tempforce command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix cac/tempforce does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix cac/tempforce command");
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;

  maxatom = 1;
  memory->create(sforce,maxatom,3,"tempforce:sforce");
}

/* ---------------------------------------------------------------------- */

FixCAC_Temp_Force::~FixCAC_Temp_Force()
{
  if (copymode) return;

  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixCAC_Temp_Force::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCAC_Temp_Force::init()
{
  // check variables
  if (!atom->CAC_flag) error->all(FLERR,"CAC fix styles require a CAC atom style");
  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR, "Variable name for fix cac/addforce does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR, "Variable for fix cac/addforce is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR, "Variable name for fix cac/addforce does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR, "Variable for fix cac/addforce is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR, "Variable name for fix cac/addforce does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR, "Variable for fix cac/addforce is invalid style");
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR, "Region ID for fix cac/addforce does not exist");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;


  // cannot use non-zero forces for a minimization since no energy is integrated
  // use fix addforce instead

  int flag = 0;
  if (update->whichflag == 2) {
    if (xstyle == EQUAL || xstyle == ATOM) flag = 1;
    if (ystyle == EQUAL || ystyle == ATOM) flag = 1;
    if (zstyle == EQUAL || zstyle == ATOM) flag = 1;
    if (xstyle == CONSTANT && xvalue != 0.0) flag = 1;
    if (ystyle == CONSTANT && yvalue != 0.0) flag = 1;
    if (zstyle == CONSTANT && zvalue != 0.0) flag = 1;
  }
  if (flag)
    error->all(FLERR, "Cannot use non-zero forces in an energy minimization");
}

/* ---------------------------------------------------------------------- */

void FixCAC_Temp_Force::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else
    error->all(FLERR, "Cannot use respa with cac/tempforce");
}

/* ---------------------------------------------------------------------- */

void FixCAC_Temp_Force::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCAC_Temp_Force::post_force(int vflag)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double ****nodal_forces = atom->nodal_forces;
  int nodes_per_element;
  int *nodes_count_list = atom->nodes_per_element_list;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  //
  int **element_scale = atom->element_scale;
  int *current_element_scale;
  int current_poly_count;
  int natoms;
  double amplitude;
  double boltz = force->boltz; 
  double nktv2p = force->nktv2p;
  double ssurfx, ssurfy, ssurfz, vtotal;
  double weigth_ratio;
  double Escale;
  double ax,ay,az,bx,by,bz,norm;
  double nx_x,nx_y,nx_z,ny_x,ny_y,ny_z,nz_x,nz_y,nz_z;
  //int timestep = update->ntimestep;
  //double dt = update->dt;
  //double time_factor;

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,3,"cac_addforce:sforce");
  }

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;
  double ****nodal_positions = atom->nodal_positions;
  double ****nodal_velocities = atom->nodal_velocities;


  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **node_types = atom->node_types;

  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++) {
	//element volume and element surface areas
	current_element_scale = element_scale[i];
	
      if (mask[i] & groupbit) {
        nodes_per_element = nodes_count_list[element_type[i]];
        for (int l = 0; l < poly_count[i]; l++) {
          for (int j = 0; j < nodes_per_element; j++) {
            if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;

			Escale = current_element_scale[0]; //elemernt scale just for element of Escale by Escale by Escale
			vtotal = Escale * Escale * Escale * unitcell_volume;
			natoms = Escale * Escale * Escale;				
			amplitude = 3 * boltz * (natoms-nodes_per_element) * target_temp / vtotal * amp_factor;
			
            foriginal[0] += nodal_forces[i][l][j][0];
            foriginal[1] += nodal_forces[i][l][j][1];
            foriginal[2] += nodal_forces[i][l][j][2];
	
			weigth_ratio=2/Escale;
			
			//Themal expansion force
			//The first term is determined by the element surface orientation and need to be set manually
			if (j == 0) {
			nodal_forces[i][l][j][0] += -15.420215 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][1] += -8.902866 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][2] += -12.590546 * amplitude * weigth_ratio;
			} else if (j == 1) {
			nodal_forces[i][l][j][0] += 15.420215 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][1] += -26.708569 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][2] += 0 * amplitude * weigth_ratio;
			} else if (j == 2) {
			nodal_forces[i][l][j][0] += 15.420215 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][1] += 8.902866 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][2] += -25.181074 * amplitude * weigth_ratio;
			} else if (j == 3) {
			nodal_forces[i][l][j][0] += -15.420215 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][1] += 26.708569 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][2] += -37.771608 * amplitude * weigth_ratio;
			} else if (j == 4) {
			nodal_forces[i][l][j][0] += -15.420215 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][1] += -8.902866 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][2] += 25.181074 * amplitude * weigth_ratio;
			} else if (j == 5) {
			nodal_forces[i][l][j][0] += 15.420215 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][1] += -26.708569 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][2] += 37.771608 * amplitude * weigth_ratio;
			} else if (j == 6) {
			nodal_forces[i][l][j][0] += 15.420215 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][1] += 8.902866 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][2] += 12.590546 * amplitude * weigth_ratio;
			} else if (j == 7) {
			nodal_forces[i][l][j][0] += -15.420215 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][1] += 26.708569 * amplitude * weigth_ratio;
			nodal_forces[i][l][j][2] += 0 * amplitude * weigth_ratio;
			}
			
          }
        }
      }
    }
  } 

}


/* ---------------------------------------------------------------------- */

void FixCAC_Temp_Force::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixCAC_Temp_Force::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,3,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixCAC_Temp_Force::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*3 * sizeof(double);
  return bytes;
}
