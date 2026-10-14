/******************************************************************************/
/*                                                                            */
/*  A Molecular Dynamics simulation of Lennard-Jones particles                */
/*                                                                            */
/*	This code is part of the course "Particle-based Simulations"              */
/*  taught at Eindhoven University of Technology.                             */
/*  No part of this code may be reproduced without permission of the author:  */
/*  Dr. Ir. E.A.J.F. Peters                                                   */
/*                                                                            */
/*  Dr. Ir. J.T. Padding:    version 1.1, 30/1/2013                           */
/*  Jeroen Hofman:           version 1.2, 28/7/2015                           */
/*  Dr. Ir. E.A.J.F. Peters: version 6.1, 17/9/2025                           */
/******************************************************************************/

/**
 * For the PBS molecular-dynamics assignment, the code needs to be extended:
 *
 * - Initialize vectors.type so particles get the proper type
 * - Implement bonds in initialise_bonds in file initialise.c
 * - Implement the bonded and non-bonded force in forces.c (Make forces type-dependent)
 * - Change the particle position initialization to account for bond lengths, angles and dihedrals
 * - Implement a Berendsen thermostat in dynamics.c
 * - Implement the needed on-the-fly data analysis
 *
 * Check every force you implement with the finite-difference force test in
 * force_test.c: set force_test > 0 in setparameters.c and the code compares
 * each analytical force against minus the numerical derivative of the
 * potential energy, and the virial against the derivative with respect to an
 * isotropic scaling of the box, before exiting.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "setparameters.h"
#include "initialise.h"
#include "nbrlist.h"
#include "forces.h"
#include "dynamics.h"
#include "memory.h"
#include "fileoutput.h"
#include "force_test.h"
#include "analysis.h"


/**
 * @brief Main MD simulation code. After initialization,
 * a velocity-Verlet scheme is executed for a specified number of time steps.
 *
 * The starting code simulates a fluid of identical Lennard-Jones particles;
 * the assignment extends it to a melt of united-atom n-pentane molecules
 * (CH3-CH2-CH2-CH2-CH3, with two dihedral angles per molecule).
 *
 * Tasks implemented in this assignment:
 *
 * 1. Add support for multiple particle types (CH3 and CH2) and handle non-bonded forces accordingly.
 * 2. Modify the code to include bonded interactions for n-pentane molecules.
 * 3. Change the particle position initialization to account for bond lengths, angles and dihedrals.
 * 4. Implement the Berendsen thermostat to maintain the system temperature (NVT ensemble).
 * 5. Perform on-the-fly analysis for the velocity distribution, the joint distribution of the
 *    two torsion angles (including the pentane effect) and the mean-square displacement.
 *
 * @return int 0 on success, non-zero on failure.
 */
int main(void)
{
    struct Vectors vectors;
    struct Parameters parameters;
    struct Nbrlist nbrlist;
    struct Analysis analysis;

    size_t step;
    double Ekin, Epot, time;


    // Step 1: Set the simulation parameters from input files
    set_parameters(&parameters);


    // Step 2: Allocate memory for particles, forces, and neighbor lists
    alloc_memory(
        &parameters,
        &vectors,
        &nbrlist
    );


    //Initialize the analysis if enabled
    if (parameters.analysis_on == 1 || parameters.velocity_analysis_on == 1 || parameters.profile_on == 1)
    {
        analysis_init(&parameters, &vectors, &analysis);
    }


    // Check if a restart is required. A restart file stores only positions,
    // velocities and forces, so the types and the molecular topology are
    // reconstructed here just as in a fresh initialisation.
    if (parameters.load_restart == 1)
    {
        load_restart(
            &parameters,
            &vectors
        );

        boundary_conditions(
            &parameters,
            &vectors
        );

        initialise_types(
            &parameters,
            &vectors
        );

        initialise_structure(
            &parameters,
            &vectors,
            &nbrlist
        );

        // Independent replicas of the same equilibrated state: redraw the
        // velocities from the Maxwell-Boltzmann distribution with another seed
        if (parameters.reseed_velocities > 0)
        {
            srand(parameters.reseed_velocities);
            initialise_velocities(
                &parameters,
                &vectors
            );
        }

        step = 0;
        time = 0.0;
    }
    else
    {
        
        initialise(
            &parameters,
            &vectors,
            &nbrlist,
            &step,
            &time
        );

        boundary_conditions(
            &parameters,
            &vectors
        );

        save_restart(
            &parameters,
            &vectors
        );
    }


    // Step 3: Build the neighbor list for non-bonded interactions
    build_nbrlist(
        &parameters,
        &vectors,
        &nbrlist
    );


    // Step 4: Calculate initial forces (non-bonded and bonded if implemented)
    Epot = calculate_forces(
        &parameters,
        &nbrlist,
        &vectors
    );


    // Output initial particle positions in PDB format
    record_trajectories_pdb(
        1,
        &parameters,
        &vectors,
        time
    );


    // CSV output header.
    // The last two columns are the temperatures calculated separately
    // for the CH3 and CH2 particle populations.
    fprintf(stdout,
        "step,time,Epot,Ekin,Etot\n"
    );

    
    // Main MD loop using velocity-Verlet integration
    while (step < parameters.num_dt_steps)
    {
        step++;
        time += parameters.dt;


        // Update velocities (half-step)
        Ekin = update_velocities_half_dt(
            &parameters,
            &nbrlist,
            &vectors
        );



        // Update positions
        update_positions(
            &parameters,
            &nbrlist,
            &vectors
        );


        // Apply boundary conditions
        boundary_conditions(
            &parameters,
            &vectors
        );


        // Rebuild neighbor list if needed
        update_nbrlist(
            &parameters,
            &vectors,
            &nbrlist
        );


        // Calculate forces for the current configuration (bonded forces if implemented)
        Epot = calculate_forces(
            &parameters,
            &nbrlist,
            &vectors
        );


        // Final velocity update (half-step)
        Ekin = update_velocities_half_dt(
            &parameters,
            &nbrlist,
            &vectors
        );

        fprintf(stdout,
                "%zu,%.15g,%.15g,%.15g,%.15g\n",
                step,
                time,
                Epot,
                Ekin,
                Epot + Ekin
            );


        // Output system state every 'num_dt_pdb' steps
        if (step % parameters.num_dt_pdb == 0)
        {
            record_trajectories_pdb(
                0,
                &parameters,
                &vectors,
                time
            );
        }


        // Save restart file every 'num_dt_restart' steps
        if (step % parameters.num_dt_restart == 0)
        {
            save_restart(
                &parameters,
                &vectors
            );
        }


        // On-the-fly analysis: (analysis.c)
        if (parameters.analysis_on == 1 || parameters.velocity_analysis_on == 1 || parameters.profile_on == 1)
        {
            analysis_update(
                &parameters,
                &vectors,
                &analysis
            );
        }
    }


    // Save final state
    save_restart(
        &parameters,
        &vectors
    );


    if (parameters.analysis_on == 1 || parameters.velocity_analysis_on == 1 || parameters.profile_on == 1)
    {
        analysis_write(
            &parameters,
            &analysis
        );

        analysis_free(
            &analysis
        );
    }


    // Step 5: Free memory and clean up
    free_memory(
        &vectors,
        &nbrlist
    );


    return 0;
}