/* all .c/.h (incl. vec3d.h) + README.md with build/run instructions; no binaries */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 6.6743015e-11

typedef struct Vec3D {
    double x, y, z;
} Vec3D;

/* Vector helper functions */

static inline Vec3D v3(double x, double y, double z)
{
    Vec3D v = {x, y, z};
    return v;
}

static inline Vec3D add(Vec3D a, Vec3D b)
{
    return v3(a.x + b.x,
              a.y + b.y,
              a.z + b.z);
}

static inline Vec3D sub(Vec3D a, Vec3D b)
{
    return v3(a.x - b.x,
              a.y - b.y,
              a.z - b.z);
}

static inline Vec3D scl(double s, Vec3D a)
{
    return v3(s * a.x,
              s * a.y,
              s * a.z);
}

static inline double dot(Vec3D a, Vec3D b)
{
    return a.x * b.x
         + a.y * b.y
         + a.z * b.z;
}

static inline Vec3D cross(Vec3D a, Vec3D b)
{
    return v3(a.y * b.z - a.z * b.y,
              a.z * b.x - a.x * b.z,
              a.x * b.y - a.y * b.x);
}

static inline double norm(Vec3D a)
{
    return sqrt(dot(a, a));
}


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