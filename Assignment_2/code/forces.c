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
    struct Vec3D fi = {0.0, 0.0, 0.0},
                 fk = {0.0, 0.0, 0.0};

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

        // Lengths of the two bonds forming the angle
        double rij_length = v3_norm(rij);
        double rkj_length = v3_norm(rkj);

        // cos(theta) = rij . rkj / (|rij| |rkj|)
        double cos_theta =
            v3_dot(rij, rkj) /
            (rij_length * rkj_length);

        // Numerical round-off can occasionally make cos(theta) very slightly
        // larger than 1 or smaller than -1. Clamp it before calling acos().
        if (cos_theta > 1.0)
            cos_theta = 1.0;
        else if (cos_theta < -1.0)
            cos_theta = -1.0;

        // Actual angle in radians
        double theta = acos(cos_theta);

        // Difference from the equilibrium angle
        double dtheta = theta - p_parameters->theta_0;

        // Harmonic angle potential:
        //
        // U_angle = 1/2 * k_theta * (theta - theta_0)^2
        Epot +=
            0.5 *
            p_parameters->k_theta *
            dtheta *
            dtheta;

        // sin(theta) appears in the derivative of acos().
        double sin_theta = sin(theta);

        // The geometry of a valid pentane molecule is never perfectly
        // collinear, but protect against division by zero from a degenerate
        // configuration.
        if (fabs(sin_theta) < 1.0e-12)
            continue;

        // dU/dtheta = k_theta * (theta - theta_0)
        //
        // The force on atom i is
        //
        // F_i =
        // k_theta*(theta-theta_0)/sin(theta)
        // *
        // [ r_kj/(|rij||rkj|)
        //   - cos(theta)*r_ij/|rij|^2 ]
        //
        // and the expression for atom k is obtained by exchanging i and k.

        double prefactor =
            p_parameters->k_theta *
            dtheta /
            sin_theta;

        struct Vec3D term_i_1 =
            v3_scl(
                1.0 / (rij_length * rkj_length),
                rkj
            );

        struct Vec3D term_i_2 =
            v3_scl(
                cos_theta / (rij_length * rij_length),
                rij
            );

        fi =
            v3_scl(
                prefactor,
                v3_sub(term_i_1, term_i_2)
            );

        struct Vec3D term_k_1 =
            v3_scl(
                1.0 / (rij_length * rkj_length),
                rij
            );

        struct Vec3D term_k_2 =
            v3_scl(
                cos_theta / (rkj_length * rkj_length),
                rkj
            );

        fk =
            v3_scl(
                prefactor,
                v3_sub(term_k_1, term_k_2)
            );

        // For a three-particle internal interaction the virial can be written
        // relative to the central atom j:
        //
        // W_angle = F_i . r_ij + F_k . r_kj
        p_vectors->press_vir_bnd +=
            v3_dot(fi, rij) +
            v3_dot(fk, rkj);

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

    struct Vec3D fi = {0.0, 0.0, 0.0};
    struct Vec3D fj = {0.0, 0.0, 0.0};
    struct Vec3D fk = {0.0, 0.0, 0.0};
    struct Vec3D fl = {0.0, 0.0, 0.0};

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

        // For compact notation use
        //
        // a = r_ij
        // b = r_kj
        // c = r_kl
        //
        // The two plane normals are
        //
        // n = a x b
        // m = b x c
        struct Vec3D a = rij;
        struct Vec3D b = rkj;
        struct Vec3D c = rkl;

        struct Vec3D n = v3_cross(a, b);
        struct Vec3D m = v3_cross(b, c);

        double N2 = v3_dot(n, n);
        double M2 = v3_dot(m, m);

        // A torsion angle is undefined if either of its planes is degenerate.
        // A correctly initialized pentane molecule should never encounter this,
        // but avoid division by zero in case of a pathological configuration.
        if (N2 < 1.0e-24 || M2 < 1.0e-24)
            continue;

        double N = sqrt(N2);
        double M = sqrt(M2);

        // q = cos(phi) = n.m / (|n||m|)
        double cos_phi =
            v3_dot(n, m) /
            (N * M);

        // Protect against tiny floating-point excursions outside [-1,1].
        if (cos_phi > 1.0)
            cos_phi = 1.0;
        else if (cos_phi < -1.0)
            cos_phi = -1.0;

        double cos2_phi = cos_phi * cos_phi;
        double cos3_phi = cos2_phi * cos_phi;

        // Ryckaert-Bellemans torsion potential:
        //
        // U(phi) =
        // c0 + c1*cos(phi)
        //    + c2*cos(phi)^2
        //    + c3*cos(phi)^3
        Epot +=
            p_parameters->c_0 +
            p_parameters->c_1 * cos_phi +
            p_parameters->c_2 * cos2_phi +
            p_parameters->c_3 * cos3_phi;

        // Differentiate the potential with respect to q = cos(phi):
        //
        // H = dU/dq
        //   = c1 + 2*c2*q + 3*c3*q^2
        double H =
            p_parameters->c_1 +
            2.0 * p_parameters->c_2 * cos_phi +
            3.0 * p_parameters->c_3 * cos2_phi;

        // Derivatives of q with respect to the two plane normals:
        //
        // A = dq/dn = m/(NM) - q*n/N^2
        //
        // B = dq/dm = n/(NM) - q*m/M^2
        struct Vec3D A =
            v3_sub(
                v3_scl(1.0 / (N * M), m),
                v3_scl(cos_phi / N2, n)
            );

        struct Vec3D B =
            v3_sub(
                v3_scl(1.0 / (N * M), n),
                v3_scl(cos_phi / M2, m)
            );

        // Using
        //
        // n = a x b
        // m = b x c
        //
        // the derivatives of q with respect to a, b and c are
        //
        // g_a = b x A
        // g_b = A x a + c x B
        // g_c = B x b
        struct Vec3D g_a =
            v3_cross(b, A);

        struct Vec3D g_b =
            v3_add(
                v3_cross(A, a),
                v3_cross(c, B)
            );

        struct Vec3D g_c =
            v3_cross(B, b);

        // Mapping the bond-vector derivatives back to the four atomic
        // coordinates gives
        //
        // grad_i(q) = g_a
        // grad_j(q) = -g_a - g_b
        // grad_k(q) = g_b + g_c
        // grad_l(q) = -g_c
        //
        // and F = -dU/dq * grad(q) = -H * grad(q).

        fi = v3_scl(-H, g_a);

        fj =
            v3_scl(
                H,
                v3_add(g_a, g_b)
            );

        fk =
            v3_scl(
                -H,
                v3_add(g_b, g_c)
            );

        fl = v3_scl(H, g_c);

        // Accumulate the torsional virial.
        //
        // Taking j as the reference position,
        //
        // r_i - r_j = a
        // r_k - r_j = b
        // r_l - r_j = b - c
        //
        // so W = F_i.a + F_k.b + F_l.(b-c).
        struct Vec3D rlj =
            v3_sub(b, c);

        p_vectors->press_vir_bnd +=
            v3_dot(fi, a) +
            v3_dot(fk, b) +
            v3_dot(fl, rlj);

        // Add the four torsion forces.
        f[i] = v3_add(f[i], fi);
        f[j] = v3_add(f[j], fj);
        f[k] = v3_add(f[k], fk);
        f[l] = v3_add(f[l], fl);
    }

    return Epot; // Return the potential energy due to dihedral-torsion interactions
}


// This function calculates non-bonded forces between particles using the neighbor list.
// The potential energy and forces are calculated using the Lennard-Jones potential.
double calculate_forces_nb(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D df;
    double r_cutsq, sr2, sr6, sr12, fr;
    struct DeltaR rij;
    struct Pair *nbr = p_nbrlist->nbr;
    const size_t num_nbrs = p_nbrlist->num_nbrs;
    struct Vec3D *f = p_vectors->f;
    struct Vec3D *r_pos = p_vectors->r;
    const struct Vec3D L = p_parameters->L;

    r_cutsq = p_parameters->r_cut * p_parameters->r_cut;
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
        rij.sq = v3_dot(rij.v, rij.v);

        // Compute forces if the distance is smaller than the cutoff distance
        if (rij.sq < r_cutsq)
        {
            double factor = nbr[k].factor; // 0 or 1, or the 1-4 scaling factor

            /// \todo Make the LJ parameters type-dependent (CH3 and CH2)

            // Get the particle types
            int type_i = p_vectors->type[i];
            int type_j = p_vectors->type[j];

            // Lorentz mixing rule:
            // sigma_ij = (sigma_i + sigma_j) / 2
            double sigma = 0.5 * (
                p_parameters->sigma[type_i] +
                p_parameters->sigma[type_j]
            );

            // Berthelot mixing rule:
            // epsilon_ij = sqrt(epsilon_i * epsilon_j)
            double epsilon = sqrt(
                p_parameters->epsilon[type_i] *
                p_parameters->epsilon[type_j]
            );

            double sigmasq = sigma * sigma;

            // The shifted LJ potential must use the parameters of this pair
            double sr2_cut = sigmasq / r_cutsq;
            double sr6_cut = sr2_cut * sr2_cut * sr2_cut;
            double sr12_cut = sr6_cut * sr6_cut;
            double Epot_cutoff = sr12_cut - sr6_cut;

            // The LJ powers are computed by squaring: (sigma/r)^6 and
            // (sigma/r)^12 from (sigma/r)^2, with no call to pow()
            sr2 = sigmasq / rij.sq;
            sr6 = sr2 * sr2 * sr2;
            sr12 = sr6 * sr6;

            // Calculate the potential energy, shifted to zero at the cutoff
            Epot += factor * 4.0 * epsilon * (sr12 - sr6 - Epot_cutoff);

            // Compute the force and apply it to both particles
            fr =
                factor *
                24.0 *
                epsilon *
                (2.0 * sr12 - sr6) /
                rij.sq;

            // Force divided by distance: multiplying by the vector rij then gives the force, with no sqrt needed

            p_vectors->press_vir_nb += fr * rij.sq;
            // virial contribution f_ij . r_ij

            df = v3_scl(fr, rij.v);

            // Update forces on particles i and j
            f[i] = v3_add(f[i], df);
            f[j] = v3_sub(f[j], df);
        }
    }

    return Epot; // Return the potential energy due to non-bonded interactions
}