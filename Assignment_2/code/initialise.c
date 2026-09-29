#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "random.h"
#include "initialise.h"

// This function initializes the particle types. The type of a particle is the
// index used to look up its type-dependent parameters, such as mass, epsilon
// and sigma (see pair_index in structs.h for the pair parameters).

void initialise_types(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        size_t position_in_molecule = i % 5;

        if (position_in_molecule == 0 || position_in_molecule == 4)
            p_vectors->type[i] = 0;   // CH3
        else
            p_vectors->type[i] = 1;   // CH2
    }
}



// This function initializes the bond connectivity between particles.
// This will be important for handling bonded interactions in the simulation.
void initialise_bond_connectivity(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    // Each n-pentane molecule contains five united-atom sites:
    //
    // CH3 - CH2 - CH2 - CH2 - CH3
    //
    // Therefore each molecule contains four covalent bonds.
    size_t num_molecules = p_parameters->num_part / 5;
    size_t num_bonds = 4 * num_molecules;

    struct Bond *bonds =
        (struct Bond *)malloc(num_bonds * sizeof(struct Bond));

    /// \todo Specify bonds between particles, i.e., bonds[i].i and bonds[i].j for bonded particle pairs.

    size_t bond_index = 0;

    // Construct the four consecutive bonds of every pentane molecule.
    //
    // For molecule 0:
    // 0-1, 1-2, 2-3, 3-4
    //
    // For molecule 1:
    // 5-6, 6-7, 7-8, 8-9
    //
    // and so on.
    for (size_t molecule = 0; molecule < num_molecules; ++molecule)
    {
        size_t first = 5 * molecule;

        for (size_t site = 0; site < 4; ++site)
        {
            bonds[bond_index].i = first + site;
            bonds[bond_index].j = first + site + 1;

            ++bond_index;
        }
    }

    p_vectors->num_bonds = num_bonds;
    p_vectors->bonds = bonds;
}


// This function derives the complete molecular structure from the list of
// bonds set up in initialise_bond_connectivity: every angle triplet i-j-k
// (two bonds sharing atom j) and every dihedral quadruplet i-j-k-l (three
// consecutive bonds), plus the 1-2, 1-3 and 1-4 partner lists that the
// neighbor list uses to exclude or scale non-bonded interactions.
//
// The partner lists all use the same compact layout: the partners of particle
// i are stored in pairs12[head12[i]] up to (excluding) pairs12[head12[i+1]].
// Such a list is built in three passes: count the partners of each particle,
// turn the counts into starting positions by a running sum, then fill the
// slots. No searching or sorting is needed.
void initialise_structure(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Nbrlist *p_nbrlist)
{
    initialise_bond_connectivity(p_parameters, p_vectors); // Initialize bonds

    struct Bond *bonds = p_vectors->bonds;
    size_t num_bonds = p_vectors->num_bonds;
    size_t num_part = p_parameters->num_part;

    // Build the 1-2 partner list from the bonds. Each bond i-j contributes two
    // entries: j is a partner of i and i is a partner of j.
    size_t *cnt = (size_t *)calloc(num_part + 1, sizeof(size_t));

    for (size_t i = 0; i < num_bonds; ++i)
    {
        ++cnt[bonds[i].i + 1];
        ++cnt[bonds[i].j + 1];
    }

    size_t *head12 = (size_t *)malloc((num_part + 1) * sizeof(size_t));

    head12[0] = 0;

    for (size_t i = 1; i <= num_part; ++i)
    {
        head12[i] = cnt[i] + head12[i - 1]; // running sum: where the partners of particle i start
        cnt[i] = head12[i];                 // reused as the next free slot while filling
    }

    size_t *pairs12 = (size_t *)malloc(cnt[num_part] * sizeof(size_t));

    for (size_t i = 0; i < num_bonds; ++i)
    {
        pairs12[cnt[bonds[i].i]++] = bonds[i].j;
        pairs12[cnt[bonds[i].j]++] = bonds[i].i;
    }


    // An angle i-j-k is a pair of bonds sharing the central atom j, so each
    // pair of 1-2 partners of an atom defines one angle centred on it.
    size_t num_angles = 0;

    for (size_t i = 1; i <= num_part; ++i)
    {
        size_t num_pairs = head12[i] - head12[i - 1];
        num_angles += (num_pairs * (num_pairs - 1)) / 2;
    }

    struct Angle *angles =
        (struct Angle *)malloc(num_angles * sizeof(struct Angle));

    size_t m = 0;

    for (size_t i = 0; i < num_part; ++i)
        for (size_t j = head12[i]; j < head12[i + 1]; ++j)
            for (size_t k = j + 1; k < head12[i + 1]; ++k)
            {
                angles[m].i = pairs12[j];
                angles[m].j = i;
                angles[m].k = pairs12[k];
                ++m;
            }


    // A dihedral i-j-k-l is built around a central bond j-k: i is any other
    // partner of j, and l any other partner of k.
    size_t num_dihdr = 0;

    for (size_t i = 0; i < num_bonds; ++i)
    {
        num_dihdr +=
            (head12[bonds[i].i + 1] - head12[bonds[i].i] - 1) *
            (head12[bonds[i].j + 1] - head12[bonds[i].j] - 1);
    }

    struct Dihedral *dihedrals =
        (struct Dihedral *)malloc(num_dihdr * sizeof(struct Dihedral));

    size_t k = 0;

    for (size_t i = 0; i < num_bonds; ++i)
    {
        struct Dihedral dihdr;

        dihdr.j = bonds[i].i;
        dihdr.k = bonds[i].j;

        for (size_t p = head12[dihdr.j]; p < head12[dihdr.j + 1]; ++p)
        {
            if (pairs12[p] == dihdr.k)
                continue;

            dihdr.i = pairs12[p];

            for (size_t q = head12[dihdr.k]; q < head12[dihdr.k + 1]; ++q)
            {
                if (pairs12[q] != dihdr.j)
                {
                    dihdr.l = pairs12[q];
                    dihedrals[k++] = dihdr;
                }
            }
        }
    }


    // 1-3 partners are the two outer atoms of each angle. The list is only
    // needed when the 1-3 non-bonded interaction is scaled (factor != 1).
    if (p_parameters->factor_13_nb != 1.0)
    {
        for (size_t i = 0; i <= num_part; ++i)
            cnt[i] = 0;

        for (size_t i = 0; i < num_angles; ++i)
        {
            ++cnt[angles[i].i + 1];
            ++cnt[angles[i].k + 1];
        }

        size_t *head13 =
            (size_t *)malloc((num_part + 1) * sizeof(size_t));

        head13[0] = 0;

        for (size_t i = 1; i <= num_part; ++i)
        {
            head13[i] = cnt[i] + head13[i - 1];
            cnt[i] = head13[i];
        }

        size_t *pairs13 =
            (size_t *)malloc(cnt[num_part] * sizeof(size_t));

        for (size_t i = 0; i < num_angles; ++i)
        {
            pairs13[cnt[angles[i].i]++] = angles[i].k;
            pairs13[cnt[angles[i].k]++] = angles[i].i;
        }

        p_nbrlist->head13 = head13;
        p_nbrlist->pairs13 = pairs13;
    }


    // 1-4 partners are the two outer atoms of each dihedral, built the same way.
    if (p_parameters->factor_14_nb != 1.0)
    {
        for (size_t i = 0; i <= num_part; ++i)
            cnt[i] = 0;

        for (size_t i = 0; i < num_dihdr; ++i)
        {
            ++cnt[dihedrals[i].i + 1];
            ++cnt[dihedrals[i].l + 1];
        }

        size_t *head14 =
            (size_t *)malloc((num_part + 1) * sizeof(size_t));

        head14[0] = 0;

        for (size_t i = 1; i <= num_part; ++i)
        {
            head14[i] = cnt[i] + head14[i - 1];
            cnt[i] = head14[i];
        }

        size_t *pairs14 =
            (size_t *)malloc(cnt[num_part] * sizeof(size_t));

        for (size_t i = 0; i < num_dihdr; ++i)
        {
            pairs14[cnt[dihedrals[i].i]++] = dihedrals[i].l;
            pairs14[cnt[dihedrals[i].l]++] = dihedrals[i].i;
        }

        p_nbrlist->head14 = head14;
        p_nbrlist->pairs14 = pairs14;
    }


    p_vectors->num_angles = num_angles;
    p_vectors->angles = angles;

    p_vectors->num_dihedrals = num_dihdr;
    p_vectors->dihedrals = dihedrals;

    // The 1-2 partner list is only kept if the neighbor list needs it to scale
    // the non-bonded interaction of bonded pairs.
    if (p_parameters->factor_12_nb != 1.0)
    {
        p_nbrlist->head12 = head12;
        p_nbrlist->pairs12 = pairs12;
    }
    else
    {
        free(head12);
        free(pairs12);
    }

    free(cnt);
}



// This function initializes the simulation by calling subroutines to initialize
// particle types, positions, velocities, and bond connectivity.
void initialise(
    struct Parameters *p_parameters,
    struct Vectors *p_vectors,
    struct Nbrlist *p_nbrlist,
    size_t *p_step,
    double *p_time)
{
    initialise_types(p_parameters, p_vectors);  // Initialize particle types
    initialise_structure(p_parameters, p_vectors, p_nbrlist);  // Initialize structure (bonds, angles, dihedrals)

    srand(SEED);  // Seed random number generator

    initialise_positions(p_parameters, p_vectors);  // Initialize particle positions
    initialise_velocities(p_parameters, p_vectors);  // Initialize particle velocities

    *p_step = 0;   // Initialize step to zero
    *p_time = 0.0; // Initialize time to zero

    return;
}



// This function initializes particle positions on a cubic lattice.
// Particles are placed in a grid with spacing based on the number of particles and the box dimensions.
void initialise_positions(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    /// \todo Change the initialization of the positions such that the particles
    /// form n-pentane molecules, with the correct bond lengths and bond angles

    // -------------------------------------------------------------------------
    // Special single-molecule initialization used for the B5 bonded-force test
    // -------------------------------------------------------------------------
    //
    // B5 requires one n-pentane molecule. For num_part = 5, construct one
    // CH3-CH2-CH2-CH2-CH3 chain with:
    //
    // bond length = r_0
    // bond angle  = theta_0
    // both torsions initially trans
    //
    // The molecule is placed near the centre of the simulation box.
    //
    // It will NOT be force-tested directly from this equilibrium structure.
    // B5 requires us to thermalise it first and then perform the finite-
    // difference tests from the resulting non-equilibrium configuration.
    if (p_parameters->num_part == 5)
    {
        double bond = p_parameters->r_0;

        // If the internal bond angle is theta_0, the forward zig-zag bond
        // direction makes an angle alpha = pi - theta_0 relative to +x.
        double alpha = M_PI - p_parameters->theta_0;

        double dx = bond * cos(alpha);
        double dy = bond * sin(alpha);

        // Build a planar all-trans zig-zag chain.
        //
        // Atom 0 ---- Atom 1
        //                  \
        //                   Atom 2 ---- Atom 3
        //                                      \
        //                                       Atom 4
        //
        // Successive bond vectors are:
        //
        // b1 = (bond, 0, 0)
        // b2 = (dx, dy, 0)
        // b3 = (bond, 0, 0)
        // b4 = (dx, dy, 0)
        //
        // This gives bond angles theta_0 and cos(phi) = -1 for both
        // dihedrals, i.e. the trans configuration.

        struct Vec3D local[5];

        local[0] = v3(0.0, 0.0, 0.0);

        local[1] =
            v3(
                bond,
                0.0,
                0.0
            );

        local[2] =
            v3(
                bond + dx,
                dy,
                0.0
            );

        local[3] =
            v3(
                2.0 * bond + dx,
                dy,
                0.0
            );

        local[4] =
            v3(
                2.0 * bond + 2.0 * dx,
                2.0 * dy,
                0.0
            );


        // Determine the centre of the molecule.
        struct Vec3D molecular_centre = {0.0, 0.0, 0.0};

        for (size_t i = 0; i < 5; ++i)
            molecular_centre =
                v3_add(molecular_centre, local[i]);

        molecular_centre =
            v3_scl(1.0 / 5.0, molecular_centre);


        // Centre of the simulation box.
        struct Vec3D box_centre =
            v3(
                0.5 * p_parameters->L.x,
                0.5 * p_parameters->L.y,
                0.5 * p_parameters->L.z
            );


        // Translate the complete molecule into the centre of the box.
        for (size_t i = 0; i < 5; ++i)
        {
            p_vectors->r[i] =
                v3_add(
                    v3_sub(local[i], molecular_centre),
                    box_centre
                );
        }

        return;
    }


    // -------------------------------------------------------------------------
    // Original multi-particle lattice initialization
    // -------------------------------------------------------------------------
    //
    // This is retained for now because B5 only concerns one molecule.
    // Packing many pentane molecules is addressed separately in task B6.

    struct Vec3D dr;  // Displacement vector for positioning particles
    struct Index3D n; // Number of grid cells along each axis
    double dl;        // Lattice spacing
    size_t ipart = 0; // Particle index; size_t, like the count it is compared
                      // against and the arrays it indexes

    // Calculate lattice spacing based on particle number and box dimensions
    dl = pow(
        p_parameters->L.x *
        p_parameters->L.y *
        p_parameters->L.z /
        ((double)p_parameters->num_part),
        1.0 / 3.0
    );

    n.i = (int)ceil(p_parameters->L.x / dl);
    n.j = (int)ceil(p_parameters->L.y / dl);
    n.k = (int)ceil(p_parameters->L.z / dl);

    dr.x = p_parameters->L.x / (double)n.i;
    dr.y = p_parameters->L.y / (double)n.j;
    dr.z = p_parameters->L.z / (double)n.k;

    ipart = 0;

    for (size_t i = 0; i < n.i; ++i)
        for (size_t j = 0; j < n.j; ++j)
            for (size_t k = 0; k < n.k; ++k, ++ipart)
            {
                if (ipart >= p_parameters->num_part)
                    break;

                p_vectors->r[ipart].x = (i + 0.5) * dr.x;
                p_vectors->r[ipart].y = (j + 0.5) * dr.y;
                p_vectors->r[ipart].z = (k + 0.5) * dr.z;
            }
}



// This function initializes the velocities of particles based on the Maxwell-Boltzmann distribution.
// The total momentum is also removed to ensure zero total momentum (important for stability).
void initialise_velocities(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    /// \todo Use the type-dependent mass, and remove the total momentum rather
    /// than the average velocity, once the particles have different masses

    struct Vec3D total_momentum = {0.0, 0.0, 0.0};
    double total_mass = 0.0;

    Vec3D *v = p_vectors->v;  // Pointer to particle velocities

    // Assign random velocities to each particle
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        // Determine whether this particle is CH3 (type 0) or CH2 (type 1)
        int type = p_vectors->type[i];

        // Use the mass corresponding to this particle type
        double m = p_parameters->mass[type];

        // Maxwell-Boltzmann velocity scale:
        // standard deviation of each velocity component = sqrt(kT / m)
        double sqrtktm = sqrt(p_parameters->kT / m);

        v[i] = v3_scl(
            sqrtktm,
            v3(gauss(), gauss(), gauss())
        );

        // Accumulate total momentum:
        // P_total = sum(m_i * v_i)
        total_momentum = v3_add(
            total_momentum,
            v3_scl(m, v[i])
        );

        // Accumulate total mass
        total_mass += m;
    }

    // Centre-of-mass velocity:
    // v_cm = P_total / M_total
    struct Vec3D v_cm =
        v3_scl(1.0 / total_mass, total_momentum);

    // Remove centre-of-mass velocity so that total momentum is zero
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        v[i] = v3_sub(v[i], v_cm);
    }
}