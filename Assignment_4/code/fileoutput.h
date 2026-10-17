#ifndef FILEOUTPUT_H_
#define FILEOUTPUT_H_

/**
 * @brief Output the particle positions to an xyz file, for visualization in
 * tools such as OVITO. Each line holds a particle name, the position, the
 * radius and the speed, in SI units. The second line of a frame follows the
 * extended-XYZ convention: the box (Lattice), the column names (Properties)
 * and the time (Time, from Vectors::time), so OVITO maps the columns itself.
 * @param reset 1: write new file or overwrite existing file, 0: append data
 * @param p_parameters used members: filename_xyz
 * @param p_vectors used members: time, r, radius, v
 */
void record_trajectories_xyz(int reset, struct Parameters * p_parameters, struct Vectors * p_vectors);

/**
 * @brief Save the state of the simulation to a binary restart file. The
 * parameters and the tangential contact histories are not stored; see
 * fileoutput.c.
 *
 * @param p_parameters used member: restart_out_filename
 * @param p_vectors the state to save
 */
void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Load a restart file written by @ref save_restart. Allocates the
 * particle arrays for the number of particles found in the file and overwrites
 * num_part.
 *
 * @param p_parameters used member: restart_in_filename; num_part is overwritten
 * @param p_vectors the state to fill
 */
void load_restart(struct Parameters *p_parameters, struct Vectors *p_vectors);

#endif /* FILEOUTPUT_H_ */
