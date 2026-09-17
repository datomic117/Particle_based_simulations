//Useful to avoid repeated importing of a .h file
#ifndef IO_H
#define IO_H

#include "vec3d.h"
#include <stdio.h>
#include <stdlib.h>

/*Declare all functions for input/output*/

int count_bodies(const char *filename);
int read_initial_conditions(const char *filename,
                            int N,
                            double *m,
                            Vec3D *r,
                            Vec3D *v);
void write_diagnostics_row(FILE *f, int step, double t, double K, double U,
                            Vec3D L, Vec3D Rcom);                    
void write_trajectory_frame(FILE *f, int N, int step, double t, const Vec3D *r);



//Necessary to make #ifndef work
#endif