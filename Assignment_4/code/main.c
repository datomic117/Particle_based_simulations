/******************************************************************************/ 
/*                                                                            */
/*  A DEM simulation with linear soft-sphere contact model                    */
/*                                                                            */
/*	This code is part of the course "Particle-based Simulations"              */
/*  taught at Eindhoven University of Technology.                             */
/*  No part of this code may be reproduced without permission of the author:  */
/*  Dr. Ir. E.A.J.F. Peters                                                   */
/*                                                                            */
/*  Dr. Ir. J.T. Padding:    version 1.1, 30/1/2013                           */
/*  Jeroen Hofman:           version 1.2, 28/7/2015                           */
/*  Dr. Ir. E.A.J.F. Peters: version 6.1, 12/10/2025                          */
/******************************************************************************/ 

/**
 * PBS Assignment 4: Silo Discharge (DEM) - Required Code Extensions
 * -----------------------------------------------------------------
 * The A4 notebook is the assignment text. Every piece of code you have to write
 * already exists as a function that compiles and is called from the right
 * place; the todo marker inside it says what to fill in (see the todo list of the
 * Doxygen documentation, html/index.html). Choose the run with the macro RUN
 * at the top of setparameters.c, or on the command line:
 *
 *   gcc -O3 -DRUN=1 *.c -o dem -lm   Part A, a test collision (RUN_COLLISION)
 *   gcc -O3 -DRUN=2 *.c -o dem -lm   Part B, a bed tilted by rotating gravity (RUN_SLOPE)
 *   gcc -O3 -DRUN=3 *.c -o dem -lm   Parts C and D, the silo (RUN_SILO)
 *
 * A. VERIFICATION, restitution (initialise.c: initialise_collision)
 *    One particle hits the bottom wall (A1), or two particles collide head-on
 *    (A2), without gravity. The table printed every step gives the kinetic
 *    energy before and after the collision, and with it the restitution
 *    coefficient.
 * B. SLOPE STABILITY (setparameters.c: update_parameters)
 *    The bed of the starting code settles, then gravity rotates slowly; find
 *    the angle at which the bed starts to move.
 * C. BUILDING THE SILO (walls.c: cylinder_wall, silo_bottom_wall;
 *    initialise.c: initialise_positions_silo; setparameters.c: update_parameters)
 *    A cylindrical wall, a bottom plate with a closable orifice, a fill that
 *    settles on the closed plate.
 * D. DISCHARGE AND THE BEVERLOO LAW (walls.c: count_in_silo)
 *    The orifice opens at t_open; the number of particles left in the silo,
 *    printed as the last column of the table, gives the discharge rate.
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
#include "walls.h"

/**
 * @brief main The main of the DEM code. After initialization,
 * a velocity-Verlet scheme is executed for a specified number of time steps.
 *
 * The starting code simulates spherical particles between a floor and a
 * ceiling; the assignment builds it out to the discharge of a cylindrical
 * silo, validated against the Beverloo law. Each time step first applies the
 * time-dependent part of the run (update_parameters: the rotation of gravity
 * in Part B, the opening of the orifice in Part D).
 *
 * The table printed to the screen has the columns step, time, coordination
 * number Z, potential and kinetic energy, their sum, and N_silo, the number of
 * particles in the silo (in Parts A and B, which have no silo, that is simply
 * every particle). Redirect it to a file: ./dem > run.dat
 *
 * @return int 0 if successful
 */
int main(void)
{
    struct Vectors vectors;
    struct Parameters parameters;
    struct Nbrlist nbrlist;
    struct Colllist colllist;
    size_t step = 0;
    double Ekin, Epot;

    set_parameters(&parameters);
    alloc_memory(&parameters, &vectors, &nbrlist, &colllist);
    if(parameters.load_restart == 1)
    {
        load_restart(&parameters, &vectors);
    }
    else
        initialise(&parameters, &vectors);

    // Put every particle inside the box before the cell list is built for the
    // first time. build_celllist turns a position straight into a cell index, so
    // a position outside the box would index outside the array of cells. During
    // the loop below boundary_conditions keeps this true; here it has to be
    // established, because an initial configuration or a restart file is free to
    // put a particle anywhere.
    boundary_conditions(&parameters, &vectors);


    build_nbrlist(&parameters, &vectors, &nbrlist);
    update_colllist(&parameters, &vectors, &nbrlist, &colllist);
    Epot = calculate_forces(&parameters, &nbrlist, &colllist, &vectors);
    record_trajectories_xyz(1, &parameters, &vectors);

    // Record the run and its geometry, so that the output can be interpreted on its own.
    printf("# run %d   box %g %g %g   R_cyl %g   z_silo %g   R_orifice %g   t_open %g   num_part %zu   dt %g\n",
           parameters.run, parameters.L.x, parameters.L.y, parameters.L.z, parameters.R_cyl,
           parameters.z_silo, parameters.R_orifice, parameters.t_open, parameters.num_part, parameters.dt);
    printf("#%7s %15s %15s %15s %15s %15s %8s\n", "Step", "Time", "Z", "Epot", "Ekin", "Etot", "N_silo");

    while (step < parameters.num_dt_steps) //start of the velocity-Verlet loop
    {

        step++;
        vectors.time += parameters.dt;
        update_parameters(&parameters, vectors.time);

        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors);
        update_positions(&parameters, &nbrlist, &vectors);

        update_tangential_displacements(&parameters, &vectors, &nbrlist, &colllist);
        boundary_conditions(&parameters, &vectors);
        update_nbrlist(&parameters, &vectors, &nbrlist);
        update_colllist(&parameters, &vectors, &nbrlist, &colllist);
        Epot = calculate_forces(&parameters, &nbrlist, &colllist, &vectors);
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors);

        // Print to the screen to monitor the progress of the simulation. Z is the
        // coordination number, the average number of contacts per particle.
        if (step%parameters.num_dt_printf ==0)
            printf("%8lu %15.7g %15.7g %15.7g %15.7g %15.7g %8zu\n", (long unsigned) step, vectors.time,
                   2.0*((double) vectors.num_contacts)/((double) parameters.num_part),
                   Epot, Ekin, Epot+Ekin, count_in_silo(&parameters, &vectors));
        if (step%parameters.num_dt_traj ==0) record_trajectories_xyz(0,&parameters,&vectors);
        if (step%parameters.num_dt_restart == 0) save_restart(&parameters,&vectors); 
    }
    save_restart(&parameters,&vectors);
    free_memory(&vectors, &nbrlist, &colllist);

    return 0;
}
