#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "random.h"

// Set the parameters of this simulation. This is the only place run settings
// live: there are no input files, so changing a run means editing the values
// below and recompiling.

void set_parameters(struct Parameters *p_parameters)
{
    // -------------------------------------------------------------------------
    // Particle types and non-bonded force field
    // -------------------------------------------------------------------------

    p_parameters->mass[0] = 1;   // A
    p_parameters->mass[1] = 1;   // B

    // Fluctuation-Dissipation constant
    p_parameters->gamma = 4.5;

    // Part B force switches: set to 0 to turn off a contribution while keeping
    // the rest of the DPD force calculation active.
    p_parameters->conservative_force_on = 1;
    p_parameters->dissipative_random_force_on = 1;

    p_parameters->a[0] = 25; // fixed a_ii: a-AA and a-BB
    p_parameters->a[1] = 28; // selectable a_ij: a-AB

    p_parameters->chain_length = 3; // chain length
    p_parameters->demix_initialisation = 1; //set to 0 for normal initialisation, 1 for demixed initialisation

    // -------------------------------------------------------------------------
    // Target temperature
    // -------------------------------------------------------------------------

    // Energy is expressed in kB K, so the numerical value equals the
    // temperature in Kelvin.
    p_parameters->kT = 1;
    p_parameters->sigma = sqrt(2.0 * p_parameters->gamma * p_parameters->kT);

    // -------------------------------------------------------------------------
    // Bonded force-field parameters for n-pentane
    // -------------------------------------------------------------------------

    // Harmonic bond potential:
    //
    // U_bond(r) = 1/2 * k_b * (r - r_0)^2

    p_parameters->r_0 = 0;          // equilibrium bond length, Angstrom
    p_parameters->k_b = 2;        // K / Angstrom^2

    // Not a finite-difference force test.
    p_parameters->force_test = 0;

    // -------------------------------------------------------------------------
    // Intramolecular non-bonded exclusions
    // -------------------------------------------------------------------------
    //
    // All intermolecular bonds are included in the nb

    p_parameters->factor_12_nb = 1.0;
    p_parameters->factor_13_nb = 1.0;
    p_parameters->factor_14_nb = 1.0;

    // -------------------------------------------------------------------------
    // Time step
    // -------------------------------------------------------------------------

    // 1 fs expressed in internal simulation time units. Where time is equal to sqrt(3)
    
    p_parameters->dt = 0.04;

    // timesteps as given by the run requirements

    p_parameters->num_dt_steps = 60000;
    // -------------------------------------------------------------------------
    // Simulation box
    // -------------------------------------------------------------------------
    //
    // Box length as given by the run requirements

    p_parameters->L =
        (struct Vec3D){
            20,
            10,
            10
        };

    // Number of particles as given by the run requirements
    p_parameters->num_part = 3* p_parameters-> L.x * p_parameters-> L.y * p_parameters-> L.z;


    // -------------------------------------------------------------------------
    // Simulation length
    // -------------------------------------------------------------------------


    // Save thermodynamic output every n steps.

    p_parameters->num_dt_output = 100;

    // -------------------------------------------------------------------------
    // Non-bonded cutoff and neighbor list
    // -------------------------------------------------------------------------

    p_parameters->r_cut = 1;

    p_parameters->r_shell = 1;

    // -------------------------------------------------------------------------
    // PDB output
    // -------------------------------------------------------------------------

    // Save one PDB frame every 1000 steps.

    p_parameters->num_dt_pdb = 1000;

    snprintf(
        p_parameters->filename_pdb,
        sizeof(p_parameters->filename_pdb),
        "data/D2_N%.0f",
        p_parameters->chain_length
    );

    p_parameters->rescale_output = 1;


    // -------------------------------------------------------------------------
    // Restart input
    // -------------------------------------------------------------------------

    //Start from saved file (For B1 and B2) from scratch: 0, from restart file: 1
    p_parameters->load_restart = 0;

    snprintf(
        p_parameters->restart_in_filename,
        sizeof(p_parameters->restart_in_filename),
        "restart_files/restart_D2_N%.0f.dat",
        p_parameters->chain_length
    );


    // -------------------------------------------------------------------------
    // Restart output
    // -------------------------------------------------------------------------

    // Save the state every 1000 steps.

    p_parameters->num_dt_restart = 1000;

    snprintf(
        p_parameters->restart_out_filename,
        sizeof(p_parameters->restart_out_filename),
        "restart_files/restart_D2_N%.0f.dat",
        p_parameters->chain_length
    );


    // -------------------------------------------------------------------------
    // On-the-fly analysis
    // -------------------------------------------------------------------------
    //
    // Velocity statistics skip the first 20% of the run and then sample every
    // 100 steps. The output prefix reflects the conservative-force setting.

    p_parameters->reseed_velocities = 0;     // 0 = keep the velocities of the restart file

    p_parameters->analysis_on = 0;
    p_parameters->velocity_analysis_on = 0;
    p_parameters->chain_size_analysis_on = 0;
    p_parameters->composition_profile_on = 1;
    p_parameters->num_dt_phi = 100;
    p_parameters->num_dt_velocity = 100;
    p_parameters->velocity_start_step = p_parameters->num_dt_steps / 5;
    p_parameters->num_dt_chain_size = 100;
    p_parameters->chain_size_start_step = p_parameters->num_dt_steps / 5;
    p_parameters->num_dt_composition = 100;
    p_parameters->num_dt_composition_block = 5000;
    p_parameters->composition_bin_width = 0.25;
    p_parameters->num_dt_msd = 10;
    p_parameters->num_dt_msd_origin = 1000;
    p_parameters->num_dt_block = 100000;

    if (p_parameters->composition_profile_on)
    {
        snprintf(
            p_parameters->filename_analysis,
            sizeof(p_parameters->filename_analysis),
            "data/D2_N%.0f",
            p_parameters->chain_length
        );
    }
    else
    {
        snprintf(
            p_parameters->filename_analysis,
            sizeof(p_parameters->filename_analysis),
            "data/B3_velocity_%s",
            p_parameters->conservative_force_on ? "on" : "off"
        );
    }


    // -------------------------------------------------------------------------
    // Minimum-image safety checks
    // -------------------------------------------------------------------------
    //
    // The cutoff must satisfy:
    //
    // r_cut <= L / 2
    //
    // so that a particle does not interact with more than one periodic image
    // of the same neighbor.

    if (p_parameters->r_cut > p_parameters->L.x / 2.0)
        fprintf(stderr, "Warning! r_cut > Lx/2");

    if (p_parameters->r_cut > p_parameters->L.y / 2.0)
        fprintf(stderr, "Warning! r_cut > Ly/2");

    if (p_parameters->r_cut > p_parameters->L.z / 2.0)
        fprintf(stderr, "Warning! r_cut > Lz/2");
}