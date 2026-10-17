#ifndef FORCES_H_
#define FORCES_H_

/**
 * @brief Calculate all forces and torques on the particles: gravity,
 * particle-particle contacts and particle-wall contacts.
 *
 * The force of each particle is initialized to gravity and its torque to zero;
 * the two contact routines then *add* their contributions and *return* their
 * potential energy. Follow this convention when adding an interaction.
 *
 * @param p_parameters used members: g, A, freq, num_part, and the contact parameters
 * @param p_nbrlist used members: nbr, tij
 * @param p_colllist used members: num_contacts, contact, and the particle-wall contact data
 * @param[out] p_vectors used members: r, v, omega, f, T
 * @return double potential energy of gravity-free conservative interactions (the contact springs)
 */
double calculate_forces(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Colllist *p_colllist, struct Vectors *p_vectors);

/**
 * @brief Calculate the particle-particle contact forces and torques with the
 * linear soft-sphere model: a spring-dashpot force in the normal direction, a
 * spring-dashpot force on the accumulated tangential displacement in the
 * tangential direction, capped by Coulomb friction when the contact slides.
 *
 * @param p_parameters used members: k_n_pp, eta_n_pp, k_t_pp, eta_t_pp, fric_pp, mass_ref
 * @param p_nbrlist used members: nbr, tij
 * @param p_colllist used members: num_contacts, contact
 * @param[out] p_vectors used members: radius, mass, v, omega, f, T
 * @return double potential energy stored in the contact springs
 */
double calculate_forces_pp(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Colllist *p_colllist, struct Vectors *p_vectors);

/**
 * @brief Calculate the particle-wall contact forces and torques, with the same
 * soft-sphere model as @ref calculate_forces_pp; the wall acts as an
 * infinitely heavy flat particle, with per-wall contact parameters.
 *
 * @param p_parameters used members: k_n_pw, eta_n_pw, k_t_pw, eta_t_pw, fric_pw, mass_ref
 * @param p_colllist used members: num_w, indcs_w, wall_id, riw, vw, tiw
 * @param[out] p_vectors used members: radius, mass, v, omega, f, T
 * @return double potential energy stored in the contact springs
 */
double calculate_forces_pw(struct Parameters *p_parameters, struct Colllist *p_colllist, struct Vectors *p_vectors);


#endif /* FORCES_H_ */
