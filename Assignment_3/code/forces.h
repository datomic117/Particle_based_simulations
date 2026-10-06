#ifndef FORCES_H_
#define FORCES_H_

#include "structs.h"

/**
 * @brief Calculate total forces on particles, including both bonded and non-bonded forces.
 *
 * This function zeroes the force array and the virial accumulators, then adds
 * the contributions of every interaction: non-bonded (Lennard-Jones) and
 * bonded (bond-stretch, angle-bend, dihedral-torsion). Each contribution
 * *adds* to the forces and virials and *returns* its potential energy; follow
 * this convention when adding an interaction.
 *
 * @param p_parameters used members: num_part, L
 * @param p_nbrlist used members: num_nbrs, nbr
 * @param p_vectors used members: r, f, press_vir_nb, press_vir_bnd
 * @return double potential energy from all interactions
 */
double calculate_forces(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Calculate non-bonded forces (e.g., Lennard-Jones interactions) between particles.
 * 
 * The non-bonded interaction is the Lennard-Jones potential, truncated at the
 * cut-off distance r_cut and shifted so the energy is zero there (an untruncated
 * potential would make the energy jump every time a pair crosses r_cut). Only
 * the pairs in the neighbor list are visited. The Lorentz-Berthelot mixing
 * rules can provide the parameters between unlike types in a multi-component
 * system.
 *
 * @param p_parameters used members: r_cut, epsilon, sigma, L
 * @param p_nbrlist used members: num_nbrs, nbr
 * @param[out] p_vectors used members: r, f, press_vir_nb
 * @return double potential energy from non-bonded interactions
 */
double calculate_forces_nb(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors);

/**
 * @brief Calculate bond-stretch forces for directly bonded particles.
 * 
 * This function calculates the bond-stretch forces between particles that are connected by
 * covalent bonds. A harmonic potential is typically used to model these interactions,
 * where the bond force is proportional to the deviation from the equilibrium bond length.
 *
 * @param p_parameters used members: r_0, k_b, L
 * @param p_vectors used members: r, f, bonds, num_bonds, press_vir_bnd
 * @return double potential energy from bond-stretch interactions
 */
double calculate_forces_bond(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Calculate angle-bend forces between sets of three bonded particles (angles).
 * 
 * This function computes the forces associated with the bending of angles between
 * bonded triplets of particles (i.e., particles i-j-k). These forces are typically modeled using
 * a harmonic potential, where the force depends on the deviation of the bond angle from its
 * equilibrium value.
 *
 * @param p_parameters used members: theta_0, k_theta, L
 * @param p_vectors used members: r, f, angles, num_angles, press_vir_bnd
 * @return double potential energy from angle-bend interactions
 */
double calculate_forces_angle(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Calculate dihedral-torsion forces between sets of four bonded particles (dihedrals).
 * 
 * Dihedral forces arise from the torsional interaction between four
 * consecutively bonded particles (i-j-k-l). The potential depends on the
 * dihedral (torsion) angle phi between the plane through i,j,k and the plane
 * through j,k,l, here expressed as a polynomial in cos(phi) (the OPLS cosine
 * series). For n-pentane this interaction sets the relative population of the
 * trans and gauche conformations.
 *
 * @param p_parameters used members: c_0 .. c_3, L
 * @param p_vectors used members: r, f, press_vir_bnd
 * @return double potential energy from dihedral-torsion interactions
 */
double calculate_forces_dihedral(struct Parameters *p_parameters, struct Vectors *p_vectors);

#endif /* FORCES_H_ */
