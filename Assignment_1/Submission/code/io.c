#include "io.h"
#include "vec3d.h"

/*File for the definition of all input/ouput functions.
This also includes the diagnostics function from B6*/

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


//Compute K, U, E, L, R_com per step and save to a new file
void write_diagnostics_row(FILE *f, int step, double t, int N, const double *m, const Vec3D *r, const Vec3D *v, double U) 
{

/* Kinetic energy */
double K = 0.0;

for (int i = 0; i < N; i++) {
    K += 0.5 * m[i] * v_dot(v[i], v[i]);
}

/* Total energy */
double E = K + U;

/* Angular momentum */
Vec3D L = v_out(0.0, 0.0, 0.0);

for (int i = 0; i < N; i++) {
    L = v_add(
        L,
        v_cross(r[i], V_scl(m[i], v[i]))
    );
}

/* Centre of mass */
Vec3D Rcom = v_out(0.0, 0.0, 0.0);
double total_mass = 0.0;

for (int i = 0; i < N; i++) {
    Rcom = v_add(
        Rcom,
        V_scl(m[i], r[i])
    );

    total_mass += m[i];
}

Rcom = V_scl(1.0 / total_mass, Rcom);


}   
                            


//Save the trajectory of the bodies every step
void write_trajectory_frame(FILE *f, int N, int step, double t, const Vec3D *r)

{


}