#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forces.h"
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"

// This function calculates all forces acting on the particles (bonded and non-bonded).
// It initializes the forces array, then calculates bond-stretch, angle-bend, dihedral-torsion,
// and non-bonded forces. The total potential energy is returned.
double calculate_forces(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D *f = p_vectors->f;
    size_t num_part = p_parameters->num_part;

    // Initialize the forces to zero for all particles
    for (size_t i = 0; i < num_part; i++)
        f[i] = (struct Vec3D){0.0, 0.0, 0.0};

    // Initialize the virial accumulators; each force routine adds its f.r contributions
    p_vectors->press_vir_bnd = 0.0;
    p_vectors->press_vir_nb = 0.0;

    // Calculate the forces and accumulate the potential energy from each type of interaction
    double Epot = 0.0;
    Epot += calculate_forces_bond(p_parameters, p_vectors);
    Epot += calculate_forces_angle(p_parameters, p_vectors);
    Epot += calculate_forces_dihedral(p_parameters, p_vectors);
    Epot += calculate_forces_nb(p_parameters, p_nbrlist, p_vectors);

    // Normalize the accumulated virials W = sum(f.r) to pressure contributions W/(3V)
    double V = p_parameters->L.x * p_parameters->L.y * p_parameters->L.z;
    p_vectors->press_vir_bnd /= (3.0 * V);
    p_vectors->press_vir_nb /= (3.0 * V);

    return Epot;
}

// This function calculates bond-stretch forces based on the current positions of the bonded particles.
// It applies the minimum image convention to calculate the distance between bonded pairs and then
// computes the force and potential energy due to the bond interaction.
double calculate_forces_bond(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    struct Bond *bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij;
    struct Vec3D fi = {0.0, 0.0, 0.0};

    // Loop through each bond and calculate the forces
    for (size_t q = 0; q < num_bonds; ++q)
    {
        size_t i = bonds[q].i;
        size_t j = bonds[q].j;

        // Apply the minimum image convention for calculating distances
        rij = v3_min_image(v3_sub(r[i], r[j]), L);

        /// \todo Provide the bond force calculation and assign forces to particles i and j
        // Newton's third law: particle j receives the opposite force
        f[i] = v3_add(f[i], fi);
        f[j] = v3_sub(f[j], fi);
    }

    return Epot; // Return the potential energy due to bond-stretch interactions
}

// This function calculates angle-bend forces based on the current positions of the angle-defined particles.
// It uses the minimum image convention and computes forces due to angle interactions.
double calculate_forces_angle(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    struct Angle *angles = p_vectors->angles;
    size_t num_angles = p_vectors->num_angles;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij, rkj;
    struct Vec3D fi = {0.0, 0.0, 0.0}, fk = {0.0, 0.0, 0.0};

    // Loop through each angle and calculate the forces
    // Note: The angles triplets ijk are computed from the bond information during initialization. This is already implemented.
    for (size_t q = 0; q < num_angles; ++q)
    {
        size_t i = angles[q].i;
        size_t j = angles[q].j;
        size_t k = angles[q].k;

        // Apply the minimum image convention for calculating distances
        rij = v3_sub(r[i], r[j]);
        rij = v3_min_image(rij, L);
        rkj = v3_sub(r[k], r[j]);
        rkj = v3_min_image(rkj, L);

        /// \todo Provide the angle force calculation and assign forces to particles i, j, and k

        // The central atom takes the opposite of both outer forces, so the
        // total force of the angle interaction is zero
        f[i] = v3_add(f[i], fi);
        f[j] = v3_sub(f[j], fi);
        f[j] = v3_sub(f[j], fk);
        f[k] = v3_add(f[k], fk);
    }
    return Epot; // Return the potential energy due to angle-bend interactions
}

// This function calculates dihedral-torsion forces based on the positions of four connected particles.
// It uses the minimum image convention and computes the forces resulting from the dihedral-torsion interaction.
double calculate_forces_dihedral(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    double Epot = 0.0;
    struct Dihedral *dihedrals = p_vectors->dihedrals;
    size_t num_dihedrals = p_vectors->num_dihedrals;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D L = p_parameters->L;
    struct Vec3D rij, rkj, rkl;
    struct Vec3D fi = {0.0, 0.0, 0.0}, fk = {0.0, 0.0, 0.0}, fl = {0.0, 0.0, 0.0};

    // Loop through each dihedral and calculate the forces
    // Note: The dihedrals are computed from the bond information during initialization. This is already implemented.
    for (size_t q = 0; q < num_dihedrals; ++q)
    {
        size_t i = dihedrals[q].i;
        size_t j = dihedrals[q].j;
        size_t k = dihedrals[q].k;
        size_t l = dihedrals[q].l;

        // Naming: r_ab = r[a] - r[b]; rij, rkj, rkl point j->i, j->k, l->k.
        // Textbook b1,b2,b3 along i->j->k->l = -rij, rkj, -rkl (signs cancel).
        // Apply the minimum image convention for calculating distances
        rij = v3_sub(r[i], r[j]);
        rij = v3_min_image(rij, L);
        rkj = v3_sub(r[k], r[j]);
        rkj = v3_min_image(rkj, L);
        rkl = v3_sub(r[k], r[l]);
        rkl = v3_min_image(rkl, L);

        /// \todo Provide the dihedral-torsion force calculation and assign forces to particles i, j, k, and l
    }

    return Epot; // Return the potential energy due to dihedral-torsion interactions
}

// This function calculates non-bonded forces between particles using the neighbor list.
// The potential energy and forces are calculated using the Lennard-Jones potential.
double calculate_forces_nb(struct Parameters *p_parameters,
                          struct Nbrlist *p_nbrlist,
                          struct Vectors *p_vectors)
{
    double Epot = 0.0;

    const double r_cutsq =
        p_parameters->r_cut * p_parameters->r_cut;

    const struct Vec3D L = p_parameters->L;
    struct Vec3D *r = p_vectors->r;
    struct Vec3D *f = p_vectors->f;
    struct Pair *nbr = p_nbrlist->nbr;

    // Only four type combinations exist. Compute their constants once per
    // force evaluation instead of repeating square roots for every pair.
    double sigma_sq_pair[NUM_TYPES][NUM_TYPES];
    double epsilon_pair[NUM_TYPES][NUM_TYPES];
    double sc6_pair[NUM_TYPES][NUM_TYPES], sc12_pair[NUM_TYPES][NUM_TYPES];
    for (int a = 0; a < NUM_TYPES; ++a)
        for (int b = 0; b < NUM_TYPES; ++b)
        {
            const double sigma = 0.5 * (p_parameters->sigma[a] + p_parameters->sigma[b]);
            sigma_sq_pair[a][b] = sigma * sigma;
            epsilon_pair[a][b] = sqrt(p_parameters->epsilon[a] * p_parameters->epsilon[b]);
            const double sc2 = sigma_sq_pair[a][b] / r_cutsq;
            sc6_pair[a][b] = sc2 * sc2 * sc2;
            sc12_pair[a][b] = sc6_pair[a][b] * sc6_pair[a][b];
        }

    for (size_t pair = 0; pair < p_nbrlist->num_nbrs; pair++)
    {
        const size_t i = nbr[pair].i;
        const size_t j = nbr[pair].j;
        const double factor = nbr[pair].factor;

        // Excluded pairs contribute neither energy nor force.
        if (factor == 0.0)
            continue;

        const struct Vec3D rij =
            v3_min_image(v3_sub(r[i], r[j]), L);

        const double r_sq = v3_dot(rij, rij);

        if (r_sq >= r_cutsq)
            continue;

        if (r_sq == 0.0)
        {
            fprintf(stderr,
                    "Error: interacting sites %zu and %zu overlap.\n",
                    i, j);
            exit(EXIT_FAILURE);
        }

        const int type_i = p_vectors->type[i];
        const int type_j = p_vectors->type[j];

        // Lorentz-Berthelot mixing rules.
        const double epsilon = epsilon_pair[type_i][type_j];
        const double sigma_sq = sigma_sq_pair[type_i][type_j];

        // Powers of sigma/r for the current separation.
        const double sr2 = sigma_sq / r_sq;
        const double sr6 = sr2 * sr2 * sr2;
        const double sr12 = sr6 * sr6;

        // The energy shift must use this pair's parameters too.
        const double sc6 = sc6_pair[type_i][type_j];
        const double sc12 = sc12_pair[type_i][type_j];

        Epot += factor * 4.0 * epsilon
              * (sr12 - sr6 - sc12 + sc6);

        // Force on i = coefficient * (r_i - r_j).
        const double coefficient =
            factor * 24.0 * epsilon
            * (2.0 * sr12 - sr6) / r_sq;

        const struct Vec3D force_ij =
            v3_scl(coefficient, rij);

        f[i] = v3_add(f[i], force_ij);
        f[j] = v3_sub(f[j], force_ij);

        // Accumulate the non-bonded virial.
        p_vectors->press_vir_nb += coefficient * r_sq;
    }

    return Epot;
}
