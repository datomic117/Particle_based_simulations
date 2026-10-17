#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "memory.h"
#include "structs.h"

void record_trajectories_xyz(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors)
/*  Write the particle positions to an xyz file, a text format that
    visualization programs such as OVITO read directly. Besides the position,
    each line carries the particle radius and speed as extra columns. The
    second line of every frame is in the extended-XYZ convention: it names the
    columns (Properties=...) and gives the box (Lattice=...) and the time, so
    that OVITO recognises the radius and the box without a column mapping and
    offers the speed as a property to colour by. All values are in SI units
    (m, m/s, s), written with enough digits that a radius keeps five
    significant figures.
    The filename (without extension) is given by p_parameters->filename_xyz.
    If reset = 1 the data is written to the file deleting data it possibly contained.
    If reset = 0 the data is appended. */
{
  FILE *fp_traj;
  // room for the configured stem plus the extension and the terminator
  char filename[sizeof p_parameters->filename_xyz + 8];

  snprintf(filename, sizeof filename, "%s%s", p_parameters->filename_xyz, ".xyz");
  if (reset == 1)
  {
    fp_traj = fopen(filename, "w");
  }
  else
  {
    fp_traj = fopen(filename, "a");
  }

  struct Vec3D L = p_parameters->L;
  fprintf(fp_traj, "%zu\n", p_parameters->num_part);
  fprintf(fp_traj, "Lattice=\"%.8g 0 0 0 %.8g 0 0 0 %.8g\" "
                   "Properties=species:S:1:pos:R:3:radius:R:1:speed:R:1 Time=%.8g\n",
          L.x, L.y, L.z, p_vectors->time);
  struct Vec3D *r = p_vectors->r;
  struct Vec3D *v = p_vectors->v;
  double * R = p_vectors->radius;
  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
      fprintf(fp_traj, "C %.8f %.8f %.8f %.8f %.6f\n", r[i].x, r[i].y, r[i].z,
          R[i], sqrt(v[i].x*v[i].x + v[i].y*v[i].y +v[i].z*v[i].z));
  }

  fclose(fp_traj);
}

void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
/* Save the state of the simulation to a binary restart file: the time, the
   number of particles, and the radius, mass, position, velocity, angular
   velocity, force and torque arrays. Parameters are NOT stored; they are set
   up again when restarting. Note that the tangential contact histories are not
   stored either, so a restarted run is not a bit-for-bit continuation. */
{
  FILE* p_file = fopen( p_parameters->restart_out_filename, "wb");
  size_t num_part = p_parameters->num_part;
  size_t sz = num_part*sizeof(struct Vec3D);
  fwrite(&p_vectors->time, sizeof(double), 1, p_file);
  fwrite(&num_part, sizeof(size_t), 1, p_file);
  fwrite(p_vectors->radius, num_part*sizeof(double), 1, p_file);
  fwrite(p_vectors->mass, num_part*sizeof(double), 1, p_file);
  fwrite(p_vectors->r, sz, 1, p_file);
  fwrite(p_vectors->v ,sz, 1, p_file);
  fwrite(p_vectors->omega, sz, 1, p_file);
  fwrite(p_vectors->f, sz, 1, p_file);
  fwrite(p_vectors->T, sz, 1, p_file);
  fclose(p_file);
}

void load_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
/* Load a restart file written by save_restart. Allocates the particle arrays
   for the number of particles found in the file and overwrites num_part. */
{
  FILE* p_file = fopen( p_parameters->restart_in_filename, "rb" );
  if (p_file == NULL)
  {
    fprintf(stderr, "cannot open restart file '%s': set load_restart = 0 in "
                    "setparameters.c, or run once to create it\n",
            p_parameters->restart_in_filename);
    exit(1);
  }
  size_t num_part;
  fread(&p_vectors->time, sizeof(double), 1, p_file);
  fread(&num_part, sizeof(size_t), 1, p_file);
  size_t sz = num_part*sizeof(struct Vec3D);
  alloc_vectors(p_vectors,num_part);
  p_parameters->num_part = num_part;
  fread(p_vectors->radius, num_part*sizeof(double), 1, p_file);
  fread(p_vectors->mass, num_part*sizeof(double), 1, p_file);
  fread(p_vectors->r, sz, 1, p_file);
  fread(p_vectors->v, sz, 1, p_file);
  fread(p_vectors->omega, sz, 1, p_file);
  fread(p_vectors->f, sz, 1, p_file);
  fread(p_vectors->T, sz, 1, p_file);
  fclose(p_file);
}
