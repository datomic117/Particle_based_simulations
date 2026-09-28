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
  // Thermal energy in the internal units derived in A1.
  p_parameters->kT = 293.0;

    // Thermostat coupling time: 0.1 ps in internal units.
  p_parameters->tau_T = 0.0911836;

  // Masses in atomic mass units.
  p_parameters->mass[TYPE_CH3] = 15.035;
  p_parameters->mass[TYPE_CH2] = 14.027;

  // Lennard-Jones energy parameters: epsilon / k_B in Kelvin.
  p_parameters->epsilon[TYPE_CH3] = 98.0;
  p_parameters->epsilon[TYPE_CH2] = 46.0;

  // Lennard-Jones length parameters in angstroms.
  p_parameters->sigma[TYPE_CH3] = 3.75;
  p_parameters->sigma[TYPE_CH2] = 3.95;
  
// The parameters below control core functionalities of the code, but many values will need to be changed
//
// Run modes (the defaults below are a production NVT run):
//  - production:      force_test = 0, num_dt_steps > 0, is_NVT = 1
//  - force check:     force_test > 0 (finite-difference test of the forces on every
//                     force_test-th particle and of the virial, after which the
//                     program exits; num_dt_steps is ignored)
//  - energy check:    is_NVT = 0 (NVE: thermostat off, so Etot must be conserved)
   p_parameters->num_part = 2000;
  p_parameters->num_dt_steps = 20000;
  p_parameters->force_test = 0;
  p_parameters->is_NVT = 1;
  p_parameters->num_dt_output = 100;
  // Scaling of the non-bonded interaction between nearby atoms of one chain:
  // 0 leaves the pair out, 1 treats it like any other pair and skips the test.
  p_parameters->factor_12_nb = 0.0;          // 1-2 connected atoms excluded from non-bonded interactions
  p_parameters->factor_13_nb = 0.0;          // 1-3 connected atoms excluded from non-bonded interactions
  p_parameters->factor_14_nb = 0.0;          // 1-4 connected atoms excluded, as TraPPE prescribes
/// \todo Set the time step, box size and cut-off distance to values appropriate
/// for n-pentane at a mass density of 626 kg/m3
   // Timestep: 1 fs in internal units.
  p_parameters->dt = 0.000911836;

  // Calculate the box size for a mass density of 626 kg/m^3.
  const double molecule_mass =
      2.0 * p_parameters->mass[TYPE_CH3]
    + 3.0 * p_parameters->mass[TYPE_CH2];

  const double total_mass_kg =
      (p_parameters->num_part / 5) * molecule_mass * 1.66053906660e-27;

  const double box_length_angstrom =
      cbrt(total_mass_kg / 626.0) / 1.0e-10;

  p_parameters->L = (struct Vec3D){
      box_length_angstrom,
      box_length_angstrom,
      box_length_angstrom
  };

  p_parameters->r_cut = 14.0;
  p_parameters->r_shell = 2.0;
  p_parameters->num_dt_pdb = 500;                           //number of time steps in between pdb outputs
  strcpy(p_parameters->filename_pdb, "../data/b3_equil");     //filename (without extension) for pdb file
  p_parameters->rescale_output = 1;                         //factor used to rescale output lengthscale (Most visualisation programs identify bonds based on distances of order 1)
  p_parameters->load_restart = 0;                           //if equal 1 restart file is loaded
  strcpy(p_parameters->restart_in_filename, "restart.dat"); //filename for loaded restart file
  p_parameters->num_dt_restart = 1000;                      // number of time steps between saves
strcpy(p_parameters->restart_out_filename, "../data/b3_equil_restart.dat"); //filename for saved restart file

  // The minimum image convention requires r_cut <= L/2: beyond that a particle
  // would interact with two images of the same neighbor.
  if (p_parameters->r_cut > p_parameters->L.x / 2.0)
    fprintf(stderr, "Warning! r_cut > Lx/2");
  if (p_parameters->r_cut > p_parameters->L.y / 2.0)
    fprintf(stderr, "Warning! r_cut > Ly/2");
  if (p_parameters->r_cut > p_parameters->L.z / 2.0)
    fprintf(stderr, "Warning! r_cut > Lz/2");
}
