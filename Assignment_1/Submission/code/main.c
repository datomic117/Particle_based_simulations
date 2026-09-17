#include "vec3d.h"
#include "io.h"
#include "forces.h"
#include "integrate.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*Argc and argv are C magic. They are basically preset outside of
the code when a command-line input is detected, with 
argc being the number of command-line tokens, and 
argv being a string of the input, which can then be interpreted
as different types dpending on the strcmp()*/
int main(int argc, char *argv[])
{
    const char *filename = "../data/bodies_2026-09-01.dat";

    /* Runtime flags: configure the run from the command line, e.g.
    "./sim --integrator euler --dt 3600 --steps 10 --trajectory-every 2"
    (defaults: velocity-Verlet, dt = 3600 s, 10 steps, save every 2nd frame) */
    Integrator_type integrator = VERLET;
    double dt = 3600.0;
    int n_steps = 10;
    int trajectory_every = 2;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--integrator") == 0 && i + 1 < argc) {
            i++;

            if (strcmp(argv[i], "euler") == 0) {
                integrator = EULER;
            } else if (strcmp(argv[i], "verlet") == 0) {
                integrator = VERLET;
            } else {
                printf("Unknown integrator: %s\n", argv[i]);
                return 1;
            }
        } else if (strcmp(argv[i], "--dt") == 0 && i + 1 < argc) {
            i++;
            //atof is for string to float
            if (atof(argv[i]) <= 0) {
                printf("Invalid --dt input!\n");
                return 1;  
            } else {
                dt = atof(argv[i]);
            }
        } else if (strcmp(argv[i], "--steps") == 0 && i + 1 < argc) {
            i++;
            //atoi is for string to integer
            if (atoi(argv[i]) <= 0) {
                printf("Invalid --steps input!\n");
                return 1;  
            } else {
                n_steps = atoi(argv[i]);
            }
        } else if (strcmp(argv[i], "--trajectory-every") == 0 && i + 1 < argc) {
            i++;
            if (atoi(argv[i]) <= 0){
                printf("Invalid --trajectory-every input!\n");
                return 1;  
            } else {
                trajectory_every = atoi(argv[i]);
            }


        }
    }

    /* Determine the number of bodies from the input file */
    int N = count_bodies(filename);

    if (N <= 0) {
        return 1;
    }

    printf("Number of bodies: %d\n", N);


    /* Dynamically allocate the simulation arrays */
    double *m = malloc(N * sizeof(double));
    Vec3D *r = malloc(N * sizeof(Vec3D));
    Vec3D *v = malloc(N * sizeof(Vec3D));
    Vec3D *a = malloc(N * sizeof(Vec3D));

    if (m == NULL || r == NULL || v == NULL || a == NULL) {
        printf("Memory allocation failed\n");

        free(m);
        free(r);
        free(v);
        free(a);

        return 1;
    }


    /* Read the initial conditions */
    if (!read_initial_conditions(filename, N, m, r, v)) {

        free(m);
        free(r);
        free(v);
        free(a);

        return 1;
    }

    printf("Initial conditions read successfully.\n");


    /* Open/create the output files (B6 diagnostics, B7 trajectory) */
    FILE *diag_file = fopen("../data/diagnostics.csv", "w");
    FILE *traj_file = fopen("../data/trajectory.csv", "w");

    if (diag_file == NULL || traj_file == NULL) {
        printf("Could not open output files\n");

        if (diag_file != NULL) fclose(diag_file);
        if (traj_file != NULL) fclose(traj_file);

        free(m);
        free(r);
        free(v);
        free(a);

        return 1;
    }

    fprintf(diag_file, "step,t,K,U,E,Lx,Ly,Lz,Rcom_x,Rcom_y,Rcom_z\n");
    fprintf(traj_file, "step,t,body,x,y,z\n");


    /* Compute the initial accelerations before stepping */
    double U = compute_acc_and_potential(N, m, r, a);

    /* Diagnostics and trajectory frame for the initial state (step 0) */
    write_diagnostics_row(diag_file, 0, 0.0, N, m, r, v, U);
    write_trajectory_frame(traj_file, N, 0, 0.0, r);

    for (int step = 1; step <= n_steps; step++) {
        U = integrate_step(integrator, N, m, r, v, a, dt);

        double t = step * dt;

        write_diagnostics_row(diag_file, step, t, N, m, r, v, U);

        if (step % trajectory_every == 0) {
            write_trajectory_frame(traj_file, N, step, t, r);
        }
    }

    printf("Ran %d steps.\n", n_steps);

    fclose(diag_file);
    fclose(traj_file);


    /* Free dynamically allocated memory */
    free(m);
    free(r);
    free(v);
    free(a);

    return 0;
}