#ifndef DYNAMICS_H_
#define DYNAMICS_H_

/**
 * @brief Update particle positions using the velocities: the drift step
 * r := r + v dt of the velocity-Verlet scheme. Also accumulates the
 * displacement of each particle since the last neighbor-list build, which
 * @ref update_nbrlist uses to decide when to rebuild.
 * @param[in] p_parameters member: dt
 * @param[out] p_nbrlist used members: dr
 * @param[in,out] p_vectors members r, dr, v
 */
void update_positions(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Update velocities by half a time step using the current forces:
 * v := v + (f/m) dt/2. Called twice per time step, before and after the
 * position update, together making up the velocity-Verlet scheme.
 * @param[in] p_parameters used members: mass, dt
 * @param[in] p_nbrlist unused
 * @param[in, out] p_vectors used members: v, f, press_kin
 * @return double kinetic energy
 */
double update_velocities_half_dt(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Apply boundary conditions: particles folded back in periodic box.
 * 
 * Ensures that particles moving beyond the simulation box
 * boundaries are wrapped around to the opposite side.
 * This maintains the periodic nature of the simulation box.
 * 
 * @param[in] p_parameters used members: L
 * @param[in, out] p_vectors used members: r
 */
void boundary_conditions(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Apply a thermostat by manipulating the particle velocities, keeping
 * the system at the set temperature (NVT ensemble). Called once per time step,
 * between the two half-kicks, when @ref Parameters.is_NVT is 1.
 * @param[in] p_parameters used members: kT, tau, dt, num_part
 * @param[in, out] p_vectors used members: v
 * @param[in] Ekin current kinetic energy
 */
void thermostat(struct Parameters *p_parameters, struct Vectors *p_vectors, double Ekin);

#endif /* DYNAMICS_H_ */
