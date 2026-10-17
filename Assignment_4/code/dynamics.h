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
 * @brief Advance the tangential displacement of every current contact by one
 * time step: the memory of how far the two surfaces have slid, which the
 * tangential contact spring acts on. Called once per time step, after the
 * position update.
 *
 * @param[in] p_parameters used member: dt
 * @param[in] p_vectors used members: v, omega
 * @param[in,out] p_nbrlist used members: nbr, tij (particle-particle histories)
 * @param[in,out] p_colllist used members: num_contacts, contact, num_w, indcs_w, riw, vw, tiw
 */
void update_tangential_displacements(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist, struct Colllist *p_colllist);

/**
 * @brief Update velocities and angular velocities by half a time step using
 * the current forces and torques: v := v + (f/m) dt/2 and
 * omega := omega + (T/I) dt/2. Called twice per time step, before and after
 * the position update, together making up the velocity-Verlet scheme.
 * @param[in] p_parameters used members: dt, num_part
 * @param[in] p_nbrlist unused
 * @param[in, out] p_vectors used members: mass, radius, v, omega, f, T
 * @return double kinetic energy, translational plus rotational
 */
double update_velocities_half_dt(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Apply boundary conditions: particles folded back in periodic box.
 * @param[in] p_parameters used members: L
 * @param[in, out] p_vectors used members: r, dr
 */
void boundary_conditions(struct Parameters *p_parameters, struct Vectors *p_vectors);



#endif /* DYNAMICS_H_ */
