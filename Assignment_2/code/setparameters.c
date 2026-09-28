#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// Set the parameters of this simulation. This is the only place run settings
// live: there are no input files, so changing a run means editing the values
// below and recompiling.
void set_parameters(struct Parameters *p_parameters)
{
/// \todo Replace these demonstration values by the parameters of the n-pentane
/// force field: type-dependent masses and LJ parameters, the bond, angle and
/// dihedral parameters, and the thermostat relaxation time
// The values below are demonstration values in reduced Lennard-Jones units
  p_parameters->kT = 1.0;                                   //thermal energy
  p_parameters->mass = 1.0;                                 //mass of a particle
  p_parameters->epsilon = 1.0;                              //LJ interaction strength
  p_parameters->sigma = 1.0;                                //LJ particle diameter

// The parameters below control core functionalities of the code, but many values will need to be changed
//
// Run modes (the defaults below are a production NVT run):
//  - production:      force_test = 0, num_dt_steps > 0, is_NVT = 1
//  - force check:     force_test > 0 (finite-difference test of the forces on every
//                     force_test-th particle and of the virial, after which the
//                     program exits; num_dt_steps is ignored)
//  - energy check:    is_NVT = 0 (NVE: thermostat off, so Etot must be conserved)
  p_parameters->num_part = 2000;             //number of particles
  p_parameters->force_test = 0;              // if > 0, test the forces on every force_test-th particle and exit
  p_parameters->is_NVT = 1;                  // if equal 1 NVT ensemble, if equal 0 NVE ensemble
  p_parameters->num_dt_output = 10;          //number of time steps between saves of output file
  // Scaling of the non-bonded interaction between nearby atoms of one chain:
  // 0 leaves the pair out, 1 treats it like any other pair and skips the test.
  p_parameters->factor_12_nb = 0.0;          // 1-2 connected atoms excluded from non-bonded interactions
  p_parameters->factor_13_nb = 0.0;          // 1-3 connected atoms excluded from non-bonded interactions
  p_parameters->factor_14_nb = 0.0;          // 1-4 connected atoms excluded, as TraPPE prescribes
/// \todo Set the time step, box size and cut-off distance to values appropriate
/// for n-pentane at a mass density of 626 kg/m3
  p_parameters->num_dt_steps = 2000;                        //number of time steps
  p_parameters->dt = 0.01;                                  //integration time step
  p_parameters->L = (struct Vec3D){14.938, 14.938, 14.938}; //box size
  p_parameters->r_cut = 2.5;                                //cut-off distance of the non-bonded interaction
  p_parameters->r_shell = 0.4;                              //shell thickness for neighbor list
  p_parameters->num_dt_pdb = 500;                           //number of time steps in between pdb outputs
  strcpy(p_parameters->filename_pdb, "trajectories");       //filename (without extension) for pdb file
  p_parameters->rescale_output = 1;                         //factor used to rescale output lengthscale (Most visualisation programs identify bonds based on distances of order 1)
  p_parameters->load_restart = 0;                           //if equal 1 restart file is loaded
  strcpy(p_parameters->restart_in_filename, "restart.dat"); //filename for loaded restart file
  p_parameters->num_dt_restart = 1000;                      // number of time steps between saves
  strcpy(p_parameters->restart_out_filename, "restart.dat");//filename for saved restart file

  // The minimum image convention requires r_cut <= L/2: beyond that a particle
  // would interact with two images of the same neighbor.
  if (p_parameters->r_cut > p_parameters->L.x / 2.0)
    fprintf(stderr, "Warning! r_cut > Lx/2");
  if (p_parameters->r_cut > p_parameters->L.y / 2.0)
    fprintf(stderr, "Warning! r_cut > Ly/2");
  if (p_parameters->r_cut > p_parameters->L.z / 2.0)
    fprintf(stderr, "Warning! r_cut > Lz/2");
}
