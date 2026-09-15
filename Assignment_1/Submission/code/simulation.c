all .c/.h (incl. vec3d.h) + README.md with build/run instructions; no binaries



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define G 6.6743015e-11


typedef struct Vec3D {
    double x, y, z;
} Vec3D;

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


int main(void)
{
    int N = 12;

    double *m = malloc(N * sizeof(double));
    Vec3D *r = malloc(N * sizeof(Vec3D));
    Vec3D *v = malloc(N * sizeof(Vec3D));
    Vec3D *a = malloc(N * sizeof(Vec3D));

    if (m == NULL || r == NULL || v == NULL || a == NULL) {
        printf("Memory allocation failed\n");
        free(m); free(r); free(v); free(a); /* free() on NULL is a no-op, safe */
        return 1;
    }

    free(m);
    free(r);
    free(v);
    free(a);

    return 0;
}