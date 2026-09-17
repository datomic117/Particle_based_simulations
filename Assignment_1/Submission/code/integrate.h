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

//Necessary to make #ifndef work
#endif
