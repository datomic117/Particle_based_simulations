//Useful to avoid repeated importing of a .h file
#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "vec3d.h"

typedef enum{
    EULER,
    VERLET
} Integrator_type;

#include <stdio.h>
#include <stdlib.h>

/* Declare the integrator functions (B5).
Both integrators leave r, v and a updated for the new time step,
and return the potential energy U at the new positions. */

double euler_step(
    int N,
    const double *m,
    Vec3D *r, Vec3D *v,
    Vec3D *a,
    double dt);

double verlet_step(
    int N,
    const double *m,
    Vec3D *r,
    Vec3D *v,
    Vec3D *a,
    double dt);

double integrate_step(
    Integrator_type integrator,
    int N, const double *m,
    Vec3D *r,
    Vec3D *v,
    Vec3D *a,
    double dt);

//Necessary to make #ifndef work
#endif
