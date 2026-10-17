#ifndef INITIALISE_H_
#define INITIALISE_H_

/**
 * @brief Set up a fresh starting configuration: seeds the random number
 * generator, calls initialise_particles, and then places the particles for the
 * chosen run: initialise_collision (Part A), initialise_positions (the bed of
 * Part B) or initialise_positions_silo (Parts C and D), followed by
 * initialise_velocities.
 *
 * @param p_parameters parameters used for initialization
 * @param p_vectors used members: type, radius, mass, r, v, omega
 * @see initialise_particles, initialise_positions, initialise_velocities
 */
void initialise(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Initialises the properties of the particles: uniformly distributed
 * radii between R_min and R_max, and the mass of a sphere of that radius.
 * @param p_parameters used members: density, R_min, R_max
 * @param p_vectors used members: type, radius, mass
 * @see initialise
 */
void initialise_particles(struct Parameters *p_parameters, struct Vectors *p_vectors);


/**
 * @brief Initialises positions on a cubic lattice, spaced so that even the
 * largest particles start without overlap.
 *
 * @param p_parameters used members: L
 * @param p_vectors used members: radius, r
 * @see initialise
 */
void initialise_positions(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Initialises velocities according to a Gaussian distribution with the
 * granular temperature Tg, removes the average velocity so the system as a
 * whole is at rest, and starts all particles without rotation.
 *
 * @param p_parameters used members: Tg
 * @param p_vectors used members: mass, v, omega
 * @see initialise
 */
void initialise_velocities(struct Parameters *p_parameters, struct Vectors *p_vectors);
/**
 * @brief Sets up the test collision of Part A: one particle moving straight
 * down onto the bottom wall (num_part = 1), or two particles approaching each
 * other head-on (num_part = 2), with speed v_impact, all of radius R_min and at
 * rest rotationally (Tasks A1, A2).
 *
 * @param p_parameters used members: num_part, v_impact, R_min, density, L
 * @param p_vectors used members: radius, mass, r, v, omega
 * @see initialise
 */
void initialise_collision(struct Parameters *p_parameters, struct Vectors *p_vectors);
/**
 * @brief Places the particles inside the silo, on a loose lattice clipped to
 * the cylinder above the closed bottom plate, from which they settle under
 * gravity into a random packing (Task C3).
 *
 * @param p_parameters used members: num_part, R_cyl, z_silo, L
 * @param p_vectors used members: radius, r
 * @see initialise
 */
void initialise_positions_silo(struct Parameters *p_parameters, struct Vectors *p_vectors);

#endif /* INITIALISE_H_ */
