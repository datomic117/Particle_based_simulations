#include "vec3d.h"
#include "io.h"
#include "forces.h"
#include "integrate.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 6.6743015e-11

/* Count the number of bodies in a snapshot file */

int count_bodies(const char *filename)
{
    FILE *file = fopen(filename, "r");

    if (file == NULL) {
        printf("Could not open input file: %s\n", filename);
        return -1;
    }

    char line[512];
    int N = 0;

    while (fgets(line, sizeof(line), file) != NULL) {

        if (line[0] == '#') {
            continue;
        }

        if (line[0] == '\n' || line[0] == '\0') {
            continue;
        }

        N++;
    }

    fclose(file);

    return N;
}


/* Read masses, positions and velocities from a snapshot file */

int read_initial_conditions(const char *filename,
                            int N,
                            double *m,
                            Vec3D *r,
                            Vec3D *v)
{
    FILE *file = fopen(filename, "r");

    if (file == NULL) {
        printf("Could not open input file: %s\n", filename);
        return 0;
    }

    char line[512];
    int i = 0;

    while (fgets(line, sizeof(line), file) != NULL) {

        if (line[0] == '#') {
            continue;
        }

        if (line[0] == '\n' || line[0] == '\0') {
            continue;
        }

        int id;
        char name[64];

        int items = sscanf(line,
                           "%d %63s %lf %lf %lf %lf %lf %lf %lf",
                           &id,
                           name,
                           &m[i],
                           &r[i].x,
                           &r[i].y,
                           &r[i].z,
                           &v[i].x,
                           &v[i].y,
                           &v[i].z);

        if (items != 9) {
            printf("Error reading input data: %s", line);
            fclose(file);
            return 0;
        }

        i++;

        if (i >= N) {
            break;
        }
    }

    fclose(file);

    if (i != N) {
        printf("Expected %d bodies, but read %d\n", N, i);
        return 0;
    }

    return 1;
}


int main(void)
{
    const char *filename = "C:\\Users\\20231118\\Particle_based_simulation\\Assignment_1\\Submission\\data\\bodies_2026-09-01.dat";

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