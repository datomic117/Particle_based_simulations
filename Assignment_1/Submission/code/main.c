#include "vec3d.h"
#include "io.h"
#include "forces.h"
#include "integrate.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int main(int argc, char *argv[])
{
    const char *filename = "../data/bodies_2026-09-01.dat";

    /* Runtime flag: pick the integrator from the command line,
    e.g. "./sim --integrator euler" (default: velocity-Verlet) */
    Integrator_type integrator = VERLET;

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


    /* Compute the initial accelerations before stepping */
    compute_acc_and_potential(N, m, r, a);

    /* Placeholder time-stepping loop; step count and dt will become
    command-line options once B6/B7 (diagnostics/trajectory output) are done */
    double dt = 3600.0;
    int n_steps = 10;

    for (int step = 0; step < n_steps; step++) {
        integrate_step(integrator, N, m, r, v, a, dt);
    }

    printf("Ran %d steps.\n", n_steps);


    /* Free dynamically allocated memory */
    free(m);
    free(r);
    free(v);
    free(a);

    return 0;
}