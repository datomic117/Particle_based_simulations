//Useful to avoid repeated importing of a .h file
#ifndef FORCES_H
#define FORCES_H

#include "vec3d.h"

#include <stdio.h>
#include <stdlib.h>

double compute_acc_and_potential(
    int N,
    const double *m,
    const Vec3D *r,
    Vec3D *a
);

//Necessary to make #ifndef work
#endif