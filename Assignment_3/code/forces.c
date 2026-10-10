#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "forces.h"
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"
#include "random.h"

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
    Epot += calculate_forces_nb(p_parameters, p_nbrlist, p_vectors);
    // removed the angle and dihedral components

    // Normalize the accumulated virials W = sum(f.r) to pressure contributions W/(3V)
    double V = p_parameters->L.x * p_parameters->L.y * p_parameters->L.z;
    p_vectors->press_vir_bnd /= (3.0 * V);
    p_vectors->press_vir_nb /= (3.0 * V);

    return Epot;
}


// This function calculates bond-stretch forces based on the current positions of the bonded particles.
// It applies the minimum image convention to calculate the distance between bonded pairs and then
// computes the force and potential energy due to the bond interaction.

//Unchanged from A2
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

        // Current bond length
        double bond_length = v3_norm(rij);

        // Difference from the equilibrium bond length
        double dr = bond_length - p_parameters->r_0;

        // Harmonic bond potential:
        //
        // U_bond = 1/2 * k_b * (r - r_0)^2
        Epot += 0.5 * p_parameters->k_b * dr * dr;

        // Force on particle i:
        //
        // F_i = -k_b * (r - r_0) * r_ij / r
        //
        // The scalar multiplying r_ij is therefore
        //
        // -k_b * (r - r_0) / r
        double force_factor =
            -p_parameters->k_b * dr / bond_length;

        fi = v3_scl(force_factor, rij);

        // Add the bonded virial contribution:
        //
        // W = F_ij . r_ij
        p_vectors->press_vir_bnd += v3_dot(fi, rij);

        // Newton's third law: particle j receives the opposite force
        f[i] = v3_add(f[i], fi);
        f[j] = v3_sub(f[j], fi);
    }

    return Epot; // Return the potential energy due to bond-stretch interactions
}

// This function calculates non-bonded forces between particles using the neighbor list.
// The potential energy and forces are calculated using the DPD repulsion (Changed from Lennard-Jones)

double calculate_forces_nb(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D r_hat,df,fc,fd,fr,ft,v_ij;
    double r_cut, r_ij, a_ij, fmag, sr2, sr6, sr12;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r_pos = p_vectors->r;
    const struct Vec3D L = p_parameters->L;

    r_cut = p_parameters->r_cut;
    double Epot = 0.0;
   

    // Loop through the neighbor list and calculate the forces for each particle pair
    for (size_t k = 0; k < num_nbrs; k++)
    {
        size_t i = nbr[k].i;
        size_t j = nbr[k].j;

        // The pair list holds only the pairs; the connecting vector is taken
        // from the current positions. A pair is far closer than half a box, so
        // the cheap minimum image applies.
        rij.v = v3_min_image(v3_sub(r_pos[i], r_pos[j]), L);
        double r_ij = v3_norm(rij.v);

        // Compute forces if the distance is smaller than the cutoff distance
        if (r_ij < r_cut)
        {
            double factor = nbr[k].factor; // 0 or 1, or the 1-4 scaling factor

            const int type_i = p_vectors->type[i];
            const int type_j = p_vectors->type[j];
            a_ij = (type_i == type_j) ? p_parameters->a[0] : p_parameters->a[1];

        
            double wC = 1.0 - (r_ij / r_cut);
            r_hat = v3_scl(1.0 / r_ij, rij.v);
            ft = (struct Vec3D){0.0, 0.0, 0.0};

            // Conservative force (disable with the switch in setparameters.c)
            if (p_parameters->conservative_force_on>0)
            {
                fmag = a_ij * wC;
                fc = v3_scl(fmag, r_hat);

                // Potential energy
                Epot += 0.5 * a_ij * r_cut * wC * wC;

                p_vectors->press_vir_nb += v3_dot(fc, rij.v);
                ft = v3_add(ft, fc);
            }

            // Dissipative + random force (disable together with the switch in setparameters.c)
            if (p_parameters->dissipative_random_force_on>0)
            {
                // Random force
                double wR = 1.0 - r_ij;
                double zetadt = gauss() / sqrt(p_parameters->dt);
                fr = v3_scl(p_parameters->sigma * wR * zetadt * r_ij, r_hat);

                // Dissipative force
                v_ij = v3_sub(p_vectors->v[i], p_vectors->v[j]);
                double Ddot = v3_dot(r_hat, v_ij);
                fd = v3_scl(-p_parameters->gamma * wR * wR * Ddot, r_hat);

                ft = v3_add(ft, v3_add(fd, fr));
            }

            // Update forces on particles i and j
            f[i] = v3_add(f[i], ft);
            f[j] = v3_sub(f[j], ft);
        }
    }

    return Epot; // Return the potential energy due to non-bonded interactions
}