//Useful to avoid repeated importing of a .h file
#ifndef VEC3D_H
#define VEC3D_H

#include <math.h>

/* File for the definition of struct Vec3D. 
Helper funcions are included because they are
static inline.*/

//Define struct using typedef
typedef struct Vec3D {
    double x;
    double y;
    double z;
} Vec3D;

//Define the vector helper functions

static inline Vec3D v_out(double x, double y, double z)
{
    Vec3D v = {x, y, z};
    return v;
}

static inline Vec3D v_add(Vec3D a, Vec3D b)
{
    return v_out(a.x + b.x,
              a.y + b.y,
              a.z + b.z);
}

static inline Vec3D v_sub(Vec3D a, Vec3D b)
{
    return v_out(a.x - b.x,
              a.y - b.y,
              a.z - b.z);
}

static inline Vec3D V_scl(double s, Vec3D a)
{
    return v_out(s * a.x,
              s * a.y,
              s * a.z);
}

static inline double v_dot(Vec3D a, Vec3D b)
{
    return a.x * b.x
         + a.y * b.y
         + a.z * b.z;
}

static inline Vec3D v_cross(Vec3D a, Vec3D b)
{
    return v_out(a.y * b.z - a.z * b.y,
              a.z * b.x - a.x * b.z,
              a.x * b.y - a.y * b.x);
}

static inline double v_norm(Vec3D a)
{
    //Removed dot for optimisation
    return sqrt(a.x * a.x
         + a.y * a.y
         + a.z * a.z);
}

//Necessary to make #ifndef work
#endif