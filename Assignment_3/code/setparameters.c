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
    /// \todo Replace these demonstration values by the parameters of the n-pentane
    /// force field: type-dependent masses and LJ parameters, the bond, angle and
    /// dihedral parameters, and the thermostat relaxation time

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
    p_parameters->a[1] = 25; // selectable a_ij: a-AB

    p_parameters->chain_length = 5; // chain length
    p_parameters->demix_initialisation = 0; //set to 0 for normal initialisation, 1 for demixed initialisation

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
    
    p_parameters->dt = 0.04*sqrt(3);


    // -------------------------------------------------------------------------
    // Simulation box
    // -------------------------------------------------------------------------
    //
    // Box length as given by the run requirements

    p_parameters->L =
        (struct Vec3D){
            10,
            10,
            10
        };

    // Number of particles as given by the run requirements
    p_parameters->num_part = 3* p_parameters-> L.x * p_parameters-> L.y * p_parameters-> L.z;


    // -------------------------------------------------------------------------
    // Simulation length
    // -------------------------------------------------------------------------

    // timesteps as given by the run requirements

    p_parameters->num_dt_steps = 0;


    // Save thermodynamic output every 10 steps.

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

    strcpy(
        p_parameters->filename_pdb,
        "c_pdb_prod_1"
    );

    p_parameters->rescale_output = 1;


    // -------------------------------------------------------------------------
    // Restart input
    // -------------------------------------------------------------------------

    // Start from a freshly generated packed B6 configuration.
    //
    // load_restart = 0 means this file is NOT actually loaded during this
    // B9 run, but keeping the path correct avoids confusion when restart mode
    // is enabled later.

    p_parameters->load_restart = 0;

    strcpy(
        p_parameters->restart_in_filename,
        "restart_files/restart_c_equil_1.dat"
    );


    // -------------------------------------------------------------------------
    // Restart output
    // -------------------------------------------------------------------------

    // Save the final B9 state after 20000 steps.

    p_parameters->num_dt_restart = 400000;

    strcpy(
        p_parameters->restart_out_filename,
        "restart_files/restart_c_prod_1.dat"
    );


    // -------------------------------------------------------------------------
    // On-the-fly analysis (C1 and C2)
    // -------------------------------------------------------------------------
    //
    // analysis_on = 1 switches on the dihedral statistics (C1); msd_on = 1 also
    // accumulates the mean-square displacement of the molecular centres of
    // mass (C2). All intervals are in time steps (1 step = 1 fs).
    //
    // Output files: <filename_analysis>_hist2d.csv, _dwell.csv, _phi_trace.csv
    // and (if msd_on) _msd.csv, written every num_dt_block steps and at the end.

    p_parameters->reseed_velocities = 0;     // 0 = keep the velocities of the restart file

    p_parameters->analysis_on = 0;
    p_parameters->msd_on = 0;

    p_parameters->num_dt_phi = 100;          // dihedrals every 0.1 ps
    p_parameters->num_dt_msd = 10;           // MSD sample every 10 fs
    p_parameters->num_dt_msd_origin = 1000;  // new MSD time origin every 1 ps
    p_parameters->num_dt_block = 100000;      // blocks of 10 ps

    strcpy(
        p_parameters->filename_analysis,
        "../data/c_prod_1"
    );


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