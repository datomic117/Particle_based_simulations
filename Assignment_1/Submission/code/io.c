#include "io.h"
#include "vec3d.h"
#include <string.h>
#include <ctype.h>

/* Case-insensitive string comparison (portable substitute for
strcasecmp/_stricmp, which aren't standard C) */
static int names_match(const char *a, const char *b)
{
    while (*a && *b) {
        if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) {
            return 0;
        }
        a++;
        b++;
    }

    return *a == '\0' && *b == '\0';
}

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


/* Read only the named bodies (case-insensitive match on the file's
"name" column), in whatever order they appear in the file */
int read_bodies_by_name(
    const char *filename,
    char **names,
    int n_names,
    double *m,
    Vec3D *r,
    Vec3D *v)
{
    FILE *file = fopen(filename, "r");

    if (file == NULL) {
        printf("Could not open input file: %s\n", filename);
        return -1;
    }

    char line[512];
    int i = 0;

    while (i < n_names && fgets(line, sizeof(line), file) != NULL) {

        if (line[0] == '#') {
            continue;
        }

        if (line[0] == '\n' || line[0] == '\0') {
            continue;
        }

        int id;
        char name[64];
        double mass;
        Vec3D pos, vel;

        int items = sscanf(line,
                           "%d %63s %lf %lf %lf %lf %lf %lf %lf",
                           &id,
                           name,
                           &mass,
                           &pos.x, &pos.y, &pos.z,
                           &vel.x, &vel.y, &vel.z);

        if (items != 9) {
            printf("Error reading input data: %s", line);
            fclose(file);
            return -1;
        }

        for (int k = 0; k < n_names; k++) {
            if (names_match(name, names[k])) {
                m[i] = mass;
                r[i] = pos;
                v[i] = vel;
                i++;
                break;
            }
        }
    }

    fclose(file);

    return i;
}


//Compute K, U, E, L, R_com per step and save to a new file
void write_diagnostics_row(
    FILE *f,
    int step, 
    double t, 
    int N, 
    const double *m, 
    const Vec3D *r, 
    const Vec3D *v, 
    double U)
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
            v_cross(r[i], v_scl(m[i], v[i]))
        );
    }

    /* Centre of mass */
    Vec3D Rcom = v_out(0.0, 0.0, 0.0);
    double total_mass = 0.0;

    for (int i = 0; i < N; i++) {
        Rcom = v_add(
            Rcom,
            v_scl(m[i], r[i])
        );

        total_mass += m[i];
    }

    Rcom = v_scl(1.0 / total_mass, Rcom);

    /* Append one CSV row (header is written once by the caller) */
    fprintf(f, "%d,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.10e\n",
            step, t, K, U, E,
            L.x, L.y, L.z,
            Rcom.x, Rcom.y, Rcom.z);
}


//Save the trajectory of the bodies every step
void write_trajectory_frame(
    FILE *f, 
    int N, 
    int step, 
    double t, 
    const Vec3D *r)

{
    /* One CSV row per body (header is written once by the caller) */
    for (int i = 0; i < N; i++) {
        fprintf(f, "%d,%.10e,%d,%.10e,%.10e,%.10e\n",
                step, t, i,
                r[i].x, r[i].y, r[i].z);
    }
}