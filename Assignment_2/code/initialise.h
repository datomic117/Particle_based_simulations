#ifndef INITIALISE_H_
#define INITIALISE_H_  

/**
 * @brief Initialise_bond_connectivity defines what particles are connected by bonds. This function is called by initialise_structure
 * 
 * @param p_parameters parameters used for initialization
 * @param[out] p_vectors member bonds is used to store the bond information
 * @see initialise_structure
 */
void initialise_bond_connectivity(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Initialise_structure creates bonds, angles, and dihedrals used for force computations.
 * 
 * This function calculates the 1-2, 1-3, and 1-4 connectivity:
 * 
 * - **1-2 (bonded)**: Refers to particles that are directly bonded, and their interaction is considered as bonded forces (e.g., bond-stretch).
 * - **1-3 (angle)**: Refers to particles connected by two bonds (i.e., i-j-k), which form an angle. These interactions are computed as angle-bending forces.
 * - **1-4 (dihedral)**: Refers to particles connected by three bonds (i.e., i-j-k-l), which form a dihedral angle. These are used to compute dihedral-torsion forces.
 * 
 * Depending on the factor_12_nb/13/14 parameters, these connected pairs are excluded
 * from (or scaled in) the non-bonded interactions, to avoid double-counting forces.
 *
 * @param p_parameters Parameters used for initialization.
 * @param[out] p_vectors Members bonds, angles, and dihedrals are used to store connectivity information.
 * @param[out] p_nbrlist Members head12, pairs12, head13, pairs13, head14, pairs14 store 1-2, 1-3, and 1-4 pairs.
 * @see initialise_bond_connectivity, initialise, is_connected_12, is_connected_13, is_connected_14
 */
void initialise_structure(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist);

/**
 * @brief Assign each particle its type. The type indexes the type-dependent
 * parameters, such as mass, epsilon and sigma.
 *
 * @param p_parameters parameters used for initialization
 * @param[out] p_vectors used members: type
 */
void initialise_types(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Set up a fresh starting configuration: calls initialise_types,
 * initialise_structure, initialise_positions and initialise_velocities, and
 * resets the step counter and simulation time.
 *
 * @param p_parameters parameters used for initialization
 * @param p_vectors used members: type, r, v, bonds, angles, dihedrals
 * @param p_nbrlist used members: head12, pairs12, head13, pairs13, head14, pairs14
 * @see initialise_types, initialise_structure, initialise_positions, initialise_velocities
 */
void initialise(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist, size_t *p_step, double *p_time);

/**
 * @brief Initialises positions on a cubic lattice
 *
 * @param p_parameters used members: L
 * @param p_vectors used members: r
 */
void initialise_positions(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Initialises velocities according to the Maxwell-Boltzmann
 * distribution at temperature kT, and removes the overall momentum so the
 * system as a whole is at rest.
 *
 * @param p_parameters used members: kT, mass
 * @param[out] p_vectors used members: v
 */
void initialise_velocities(struct Parameters *p_parameters, struct Vectors *p_vectors);

#endif /* INITIALISE_H_ */
