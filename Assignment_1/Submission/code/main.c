#include "vec3d.h"
#include "io.h"
#include "forces.h"
#include "integrate.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 6.6743015e-11

/* Count the number of bodies in a snapshot file */




int main(void)
{
    const char *filename = "../data/bodies_2026-09-01.dat";

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


    /* Free dynamically allocated memory */
    free(m);
    free(r);
    free(v);
    free(a);

    return 0;
}