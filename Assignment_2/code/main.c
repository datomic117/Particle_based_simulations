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
    size_t step; 
    double Ekin, Epot, time; 

    // Step 1: Set the simulation parameters from input files
    set_parameters(&parameters); 

    // Step 2: Allocate memory for particles, forces, and neighbor lists
    alloc_memory(&parameters, &vectors, &nbrlist); 

    // Check if a restart is required. A restart file stores only positions,
    // velocities and forces, so the types and the molecular topology are
    // reconstructed here just as in a fresh initialisation.
    if (parameters.load_restart == 1)
    {
        load_restart(&parameters, &vectors);
        boundary_conditions(&parameters, &vectors);
        initialise_types(&parameters, &vectors);
        initialise_structure(&parameters, &vectors, &nbrlist);
        step = 0;
        time = 0.0;
    }
    else 
    {   
    /// \todo Initialize particle types (CH3 and CH2) in vectors.type array
    /// \todo Implement the bonds between the UA of n-pentane in initialise_bonds (initialise.c)
        initialise(&parameters, &vectors, &nbrlist, &step, &time); 
        boundary_conditions(&parameters, &vectors);
    }

    // Step 3: Build the neighbor list for non-bonded interactions
    build_nbrlist(&parameters, &vectors, &nbrlist); 

    // Step 4: Calculate initial forces (non-bonded and bonded if implemented)
    Epot = calculate_forces(&parameters, &nbrlist, &vectors); 

    // Output initial particle positions in PDB format
    record_trajectories_pdb(1, &parameters, &vectors, time); 

    // In force-test mode the forces and the virial are verified against finite
    // differences of the potential energy, after which the program exits. Every
    // force_test-th particle is tested, so force_test = 1 tests all of them.
    // Use this to check every force term you implement: the analytical force
    // must equal minus the numerical derivative of the potential energy.
    //
    // forces_test takes the routine to check as its first argument, so a term
    // can be tested on its own instead of inside the sum of all of them. That
    // is what task B5 asks for, and it matters: an error in one term can hide
    // behind a larger correct one in the total. Write a wrapper per term with
    // the signature of calculate_forces, for example
    //
    //     double forces_bond_only(struct Parameters *p_parameters,
    //                             struct Nbrlist *p_nbrlist,
    //                             struct Vectors *p_vectors)
    //     {
    //         (void)p_nbrlist;                       // this term needs no pairs
    //         for (size_t i = 0; i < p_parameters->num_part; i++)
    //             p_vectors->f[i] = (struct Vec3D){0.0, 0.0, 0.0};
    //         return calculate_forces_bond(p_parameters, p_vectors);
    //     }
    //
    // and pass it here in place of calculate_forces.
    if (parameters.force_test > 0)
    {
        for(size_t i=0; i<parameters.num_part; i+=(size_t)parameters.force_test)
            forces_test(calculate_forces, (int)i, &parameters, &nbrlist, &vectors);
        virial_test(calculate_forces, &parameters, &nbrlist, &vectors);
        free_memory(&vectors, &nbrlist);
        return 0;
    }



    printf("step,time_internal,Epot_internal,Ekin_internal,Etot_internal,T_internal\n");

    // Main MD loop using velocity-Verlet integration
    while (step < parameters.num_dt_steps) 
    { 
        step++;
        time += parameters.dt; 

        // Update velocities (half-step)
        /// \todo Implement the use of type-dependent masses
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors); 

        /// \todo Implement and apply the Berendsen thermostat to maintain temperature (dynamics.c)
        if (parameters.is_NVT == 1)
            thermostat(&parameters, &vectors, Ekin); 

        // Update positions
        update_positions(&parameters, &nbrlist, &vectors); 

        // Apply boundary conditions
        boundary_conditions(&parameters, &vectors); 

        // Rebuild neighbor list if needed
        update_nbrlist(&parameters, &vectors, &nbrlist); 

        // Calculate forces for the current configuration (bonded forces if implemented)
        Epot = calculate_forces(&parameters, &nbrlist, &vectors); 

        // Final velocity update (half-step)
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors); 

        // Output system state every 'num_dt_pdb' steps
        if (step % parameters.num_dt_pdb == 0) 
            record_trajectories_pdb(0, &parameters, &vectors, time); 

        // Save restart file every 'num_dt_restart' steps
        if (step % parameters.num_dt_restart == 0) 
            save_restart(&parameters, &vectors); 

        /// \todo Implement on-the-fly analysis of velocity distribution, torsion angle distribution and mean-square displacement
        // Print to the screen to monitor the progress of the simulation
       
        if (step % parameters.num_dt_output == 0)
        {
            const double dof = 3.0 * (double)parameters.num_part - 3.0;
            const double temperature = 2.0 * Ekin / dof;

            printf("%zu,%.15g,%.15g,%.15g,%.15g,%.15g\n",
                   step, time, Epot, Ekin, Epot + Ekin, temperature);
        }




    }

    // Save final state
    save_restart(&parameters, &vectors); 

    // Step 5: Free memory and clean up
    free_memory(&vectors, &nbrlist); 

    return 0; 
}
