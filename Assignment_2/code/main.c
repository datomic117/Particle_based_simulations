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


/*
 * B5 helper routines
 *
 * forces_test() expects a function with the same arguments as
 * calculate_forces(). These wrappers allow each bonded interaction to be
 * tested separately.
 *
 * The force array must be reset before each evaluation because the individual
 * bonded routines add their contributions to the existing forces.
 */


/* Test the harmonic bond force only. */
double forces_bond_only(
    struct Parameters *p_parameters,
    struct Nbrlist *p_nbrlist,
    struct Vectors *p_vectors)
{
    (void)p_nbrlist;

    for (size_t i = 0; i < p_parameters->num_part; ++i)
        p_vectors->f[i] = (struct Vec3D){0.0, 0.0, 0.0};

    p_vectors->press_vir_bnd = 0.0;
    p_vectors->press_vir_nb = 0.0;

    return calculate_forces_bond(
        p_parameters,
        p_vectors
    );
}


/* Test the harmonic angle force only. */
double forces_angle_only(
    struct Parameters *p_parameters,
    struct Nbrlist *p_nbrlist,
    struct Vectors *p_vectors)
{
    (void)p_nbrlist;

    for (size_t i = 0; i < p_parameters->num_part; ++i)
        p_vectors->f[i] = (struct Vec3D){0.0, 0.0, 0.0};

    p_vectors->press_vir_bnd = 0.0;
    p_vectors->press_vir_nb = 0.0;

    return calculate_forces_angle(
        p_parameters,
        p_vectors
    );
}


/* Test the Ryckaert-Bellemans torsion force only. */
double forces_dihedral_only(
    struct Parameters *p_parameters,
    struct Nbrlist *p_nbrlist,
    struct Vectors *p_vectors)
{
    (void)p_nbrlist;

    for (size_t i = 0; i < p_parameters->num_part; ++i)
        p_vectors->f[i] = (struct Vec3D){0.0, 0.0, 0.0};

    p_vectors->press_vir_bnd = 0.0;
    p_vectors->press_vir_nb = 0.0;

    return calculate_forces_dihedral(
        p_parameters,
        p_vectors
    );
}


/*
 * B5 torsion conservation check.
 *
 * Each torsion is evaluated separately. For one isolated i-j-k-l torsion we
 * calculate:
 *
 *     F_net = F_i + F_j + F_k + F_l
 *
 * and
 *
 *     tau_net = sum_s x_s x F_s
 *
 * using atom j as the origin.
 *
 * The unwrapped relative coordinates are
 *
 *     x_j = 0
 *     x_i = r_ij
 *     x_k = r_kj
 *     x_l = r_kj - r_kl
 *
 * Both F_net and tau_net should vanish to floating-point round-off.
 */
void check_torsion_force_and_torque(
    struct Parameters *p_parameters,
    struct Nbrlist *p_nbrlist,
    struct Vectors *p_vectors)
{
    struct Dihedral *all_dihedrals =
        p_vectors->dihedrals;

    size_t all_num_dihedrals =
        p_vectors->num_dihedrals;

    struct Vec3D *r =
        p_vectors->r;

    struct Vec3D *f =
        p_vectors->f;

    struct Vec3D L =
        p_parameters->L;


    printf("\n");
    printf("============================================================\n");
    printf("B5 TORSION NET-FORCE AND NET-TORQUE CHECK\n");
    printf("============================================================\n");


    for (size_t q = 0; q < all_num_dihedrals; ++q)
    {
        /*
         * Temporarily make the torsion routine see only this one dihedral.
         *
         * This is useful because we want the four forces generated by one
         * individual torsion rather than the sum of both pentane torsions.
         */
        p_vectors->dihedrals =
            &all_dihedrals[q];

        p_vectors->num_dihedrals =
            1;


        /*
         * Call the actual torsion implementation.
         *
         * forces_dihedral_only() also clears the force array first.
         */
        forces_dihedral_only(
            p_parameters,
            p_nbrlist,
            p_vectors
        );


        size_t i =
            all_dihedrals[q].i;

        size_t j =
            all_dihedrals[q].j;

        size_t k =
            all_dihedrals[q].k;

        size_t l =
            all_dihedrals[q].l;


        /*
         * Reconstruct a consistent unwrapped geometry relative to atom j.
         */
        struct Vec3D rij =
            v3_min_image(
                v3_sub(r[i], r[j]),
                L
            );

        struct Vec3D rkj =
            v3_min_image(
                v3_sub(r[k], r[j]),
                L
            );

        struct Vec3D rkl =
            v3_min_image(
                v3_sub(r[k], r[l]),
                L
            );


        /*
         * Coordinates relative to atom j.
         */
        struct Vec3D xi =
            rij;

        struct Vec3D xj =
            v3(0.0, 0.0, 0.0);

        struct Vec3D xk =
            rkj;

        struct Vec3D xl =
            v3_sub(
                rkj,
                rkl
            );


        /*
         * Sum the four torsional forces.
         */
        struct Vec3D F_net =
            v3_add(
                v3_add(
                    f[i],
                    f[j]
                ),
                v3_add(
                    f[k],
                    f[l]
                )
            );


        /*
         * Net torque:
         *
         * tau = x_i x F_i
         *     + x_j x F_j
         *     + x_k x F_k
         *     + x_l x F_l
         *
         * x_j = 0, but keep the term explicitly so the definition is clear.
         */
        struct Vec3D tau_i =
            v3_cross(
                xi,
                f[i]
            );

        struct Vec3D tau_j =
            v3_cross(
                xj,
                f[j]
            );

        struct Vec3D tau_k =
            v3_cross(
                xk,
                f[k]
            );

        struct Vec3D tau_l =
            v3_cross(
                xl,
                f[l]
            );

        struct Vec3D tau_net =
            v3_add(
                v3_add(
                    tau_i,
                    tau_j
                ),
                v3_add(
                    tau_k,
                    tau_l
                )
            );


        /*
         * Relative residuals are useful because they show how small the
         * cancellation error is compared with the individual forces and
         * torques.
         */
        double force_scale =
            v3_norm(f[i]) +
            v3_norm(f[j]) +
            v3_norm(f[k]) +
            v3_norm(f[l]);

        double torque_scale =
            v3_norm(xi) * v3_norm(f[i]) +
            v3_norm(xj) * v3_norm(f[j]) +
            v3_norm(xk) * v3_norm(f[k]) +
            v3_norm(xl) * v3_norm(f[l]);


        double relative_force_residual =
            0.0;

        double relative_torque_residual =
            0.0;


        if (force_scale > 0.0)
        {
            relative_force_residual =
                v3_norm(F_net) /
                force_scale;
        }


        if (torque_scale > 0.0)
        {
            relative_torque_residual =
                v3_norm(tau_net) /
                torque_scale;
        }


        printf(
            "\nTorsion %zu: atoms %zu-%zu-%zu-%zu\n",
            q + 1,
            i,
            j,
            k,
            l
        );

        printf(
            "Sum of four forces = "
            "(%.15e, %.15e, %.15e)\n",
            F_net.x,
            F_net.y,
            F_net.z
        );

        printf(
            "|Sum of four forces| = %.15e\n",
            v3_norm(F_net)
        );

        printf(
            "Net torque = "
            "(%.15e, %.15e, %.15e)\n",
            tau_net.x,
            tau_net.y,
            tau_net.z
        );

        printf(
            "|Net torque| = %.15e\n",
            v3_norm(tau_net)
        );

        printf(
            "Relative force residual = %.15e\n",
            relative_force_residual
        );

        printf(
            "Relative torque residual = %.15e\n",
            relative_torque_residual
        );
    }


    /*
     * IMPORTANT:
     *
     * Restore the original pointer and count. free_memory() must receive the
     * original allocated dihedral array rather than a pointer into its middle.
     */
    p_vectors->dihedrals =
        all_dihedrals;

    p_vectors->num_dihedrals =
        all_num_dihedrals;
}


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
    alloc_memory(
        &parameters,
        &vectors,
        &nbrlist
    );


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

        step = 0;
        time = 0.0;
    }
    else
    {
        /// \todo Initialize particle types (CH3 and CH2) in vectors.type array
        /// \todo Implement the bonds between the UA of n-pentane in initialise_bonds (initialise.c)

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
        /*
         * B5 requires finite-difference verification separately for:
         *
         * 1. bond forces
         * 2. angle forces
         * 3. torsion forces
         *
         * All tests below use exactly the same thermalised configuration.
         */

        printf("\n");
        printf("============================================================\n");
        printf("B5 FINITE-DIFFERENCE TEST: BOND FORCE\n");
        printf("============================================================\n");

        for (
            size_t i = 0;
            i < parameters.num_part;
            i += (size_t)parameters.force_test
        )
        {
            forces_test(
                forces_bond_only,
                (int)i,
                &parameters,
                &nbrlist,
                &vectors
            );
        }


        printf("\n");
        printf("============================================================\n");
        printf("B5 FINITE-DIFFERENCE TEST: ANGLE FORCE\n");
        printf("============================================================\n");

        for (
            size_t i = 0;
            i < parameters.num_part;
            i += (size_t)parameters.force_test
        )
        {
            forces_test(
                forces_angle_only,
                (int)i,
                &parameters,
                &nbrlist,
                &vectors
            );
        }


        printf("\n");
        printf("============================================================\n");
        printf("B5 FINITE-DIFFERENCE TEST: TORSION FORCE\n");
        printf("============================================================\n");

        for (
            size_t i = 0;
            i < parameters.num_part;
            i += (size_t)parameters.force_test
        )
        {
            forces_test(
                forces_dihedral_only,
                (int)i,
                &parameters,
                &nbrlist,
                &vectors
            );
        }


        /*
         * The assignment also requires the numerical sum of the four forces
         * and their net torque for each torsion.
         */
        check_torsion_force_and_torque(
            &parameters,
            &nbrlist,
            &vectors
        );


        free_memory(
            &vectors,
            &nbrlist
        );

        return 0;
    }


    // CSV output header.
    // The last two columns are the temperatures calculated separately
    // for the CH3 and CH2 particle populations.
    printf(
        "step,time,Epot,Ekin,Etot,T,T_CH3_K,T_CH2_K\n"
    );


    // Main MD loop using velocity-Verlet integration
    while (step < parameters.num_dt_steps)
    {
        step++;
        time += parameters.dt;


        // Update velocities (half-step)
        /// \todo Implement the use of type-dependent masses
        Ekin = update_velocities_half_dt(
            &parameters,
            &nbrlist,
            &vectors
        );


        /// \todo Implement and apply the Berendsen thermostat to maintain temperature (dynamics.c)


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


        // Apply the thermostat after the complete velocity-Verlet step.
        if (parameters.is_NVT == 1)
        {
            thermostat(
                &parameters,
                &vectors,
                Ekin
            );


            // Recalculate kinetic energy after scaling
            // the velocities with the thermostat.
            Ekin = 0.0;

            for (size_t i = 0; i < parameters.num_part; i++)
            {
                const int type =
                    vectors.type[i];

                const double mass =
                    parameters.mass[type];

                Ekin +=
                    0.5
                    * mass
                    * v3_dot(
                        vectors.v[i],
                        vectors.v[i]
                    );
            }
        }


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


        /// \todo Implement on-the-fly analysis of velocity distribution, torsion angle distribution and mean-square displacement


        // Print to the screen to monitor the progress of the simulation
        if (step % parameters.num_dt_output == 0)
        {
            // Total system temperature.
            // Three degrees of freedom are removed because
            // the centre-of-mass momentum is zero.
            const double dof =
                3.0 * (double)parameters.num_part - 3.0;

            const double temperature =
                2.0 * Ekin / dof;


            // Calculate the temperature implied by each
            // particle type separately.
            //
            // type 0 = CH3
            // type 1 = CH2
            double sum_mv2[NUM_TYPES] =
                {0.0};

            size_t count[NUM_TYPES] =
                {0};


            for (size_t i = 0; i < parameters.num_part; i++)
            {
                const int type =
                    vectors.type[i];

                const double mass =
                    parameters.mass[type];

                sum_mv2[type] +=
                    mass
                    * v3_dot(
                        vectors.v[i],
                        vectors.v[i]
                    );

                count[type]++;
            }


            // For one particle type:
            //
            //     T = sum(m v^2) / (3 N)
            //
            // Energy is expressed in kB K, so the numerical
            // result is directly the temperature in Kelvin.
            const double T_CH3 =
                sum_mv2[0]
                / (3.0 * (double)count[0]);

            const double T_CH2 =
                sum_mv2[1]
                / (3.0 * (double)count[1]);


            printf(
                "%zu,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g\n",
                step,
                time,
                Epot,
                Ekin,
                Epot + Ekin,
                temperature,
                T_CH3,
                T_CH2
            );
        }
    }


    // Save final state
    save_restart(
        &parameters,
        &vectors
    );


    // Step 5: Free memory and clean up
    free_memory(
        &vectors,
        &nbrlist
    );


    return 0;
}