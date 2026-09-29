#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "memory.h"
#include "structs.h"

// Write the particle positions to a pdb file, a text format that visualization
// programs such as VMD and OVITO read directly.
// The filename (without extension) is given by p_parameters->filename_pdb.
// If reset = 1 the data is written to the file deleting data it possibly contained.
// If reset = 0 the data is appended.
void record_trajectories_pdb(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time)
{
  FILE *fp_traj;
  char filename[1029];
  double rs = p_parameters->rescale_output;

  snprintf(filename, sizeof(filename), "%s%s", p_parameters->filename_pdb, ".pdb");
  if (reset == 1)
  {
    fp_traj = fopen(filename, "w");
  }
  else
  {
    fp_traj = fopen(filename, "a");
  }

  fprintf(fp_traj, "MODEL\n");
  fprintf(fp_traj, "REMARK TIME = %f\n", time);
  fprintf(fp_traj, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-10s%-3s\n", rs*p_parameters->L.x, rs*p_parameters->L.y, rs*p_parameters->L.z, 90.0, 90.0, 90.0, "P 1", "1");
  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
    fprintf(fp_traj, "HETATM%5u  C   UNK A   1    %8.3f%8.3f%8.3f  1.00  0.00           C\n", (unsigned int)i % 100000, rs*p_vectors->r[i].x, rs*p_vectors->r[i].y, rs*p_vectors->r[i].z);
  }
  fprintf(fp_traj, "ENDMDL\n");

  fclose(fp_traj);
}

// Write the particle positions to a xyz file
// The filename (without extension) is given by p_parameters->filename_xyz.
// If reset = 1 the data is written to the file deleting data it possibly contained.
// If reset = 0 the data is appended.
void record_trajectories_xyz(int reset, struct Parameters *p_parameters, struct Vectors *p_vectors, double time)
{
  FILE *fp_traj;
  char filename[1029];
  double rs = p_parameters->rescale_output;

  snprintf(filename, sizeof(filename), "%s%s", p_parameters->filename_xyz, ".xyz");
  if (reset == 1)
  {
    fp_traj = fopen(filename, "w");
  }
  else
  {
    fp_traj = fopen(filename, "a");
  }

  fprintf(fp_traj, "%zu\n", p_parameters->num_part);
  fprintf(fp_traj, "time = %f\n", time);
  struct Vec3D *r = p_vectors->r;
  for (size_t i = 0; i < p_parameters->num_part; i++)
  {
    fprintf(fp_traj, "  C        %10.5f %10.5f %10.5f\n", rs*r[i].x, rs*r[i].y, rs*r[i].z);
  }

  fclose(fp_traj);
}

// Save the state of the simulation to a binary restart file: the number of
// particles followed by the position, velocity and force arrays. Parameters,
// types and topology are NOT stored; they are set up again when restarting
// (see main.c), so the restart file stays valid when settings change.
void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
  FILE* p_file = fopen( p_parameters->restart_out_filename, "wb");
  size_t num_part = p_parameters->num_part;
  size_t sz = num_part*sizeof(struct Vec3D);

  fwrite(&num_part, sizeof(size_t), 1, p_file);
  fwrite(p_vectors->r, sz, 1, p_file);
  fwrite(p_vectors->v ,sz, 1, p_file);
  fwrite(p_vectors->f, sz, 1, p_file);
  fclose(p_file);
}

// Load into the arrays already allocated by main. The saved particle count
// must match the configured system, including its box and neighbour storage.
void load_restart(struct Parameters *p_parameters, struct Vectors *p_vectors)
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
  if (fread(&num_part, sizeof(size_t), 1, p_file) != 1 || num_part != p_parameters->num_part)
  {
    fprintf(stderr, "Invalid restart particle count or incompatible configuration.\n");
    fclose(p_file);
    exit(EXIT_FAILURE);
  }
  size_t sz = num_part*sizeof(struct Vec3D);
  if (fread(p_vectors->r, sz, 1, p_file) != 1 ||
      fread(p_vectors->v, sz, 1, p_file) != 1 ||
      fread(p_vectors->f, sz, 1, p_file) != 1)
  {
    fprintf(stderr, "Truncated restart file.\n");
    fclose(p_file);
    exit(EXIT_FAILURE);
  }
  fclose(p_file);
}
