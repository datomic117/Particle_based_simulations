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
    size_t c_l = p_parameters->chain_length;
    size_t num_mol = p_parameters->num_part / c_l;
    for (int m=0; m<num_mol; m++)
    {
        if (m % 2 == 0)
        {
            for  (int i=0; i<c_l; i++)
            {
                p_vectors->type[m*c_l + i] = 0; // A
            }
            
        }
        else
        {
            for (int i=0; i<c_l; i++)
            {
                p_vectors->type[m*c_l + i] = 1;   // B
            }
        }     
         
    }
}



// This function initializes the bond connectivity between particles.
// This will be important for handling bonded interactions in the simulation.
void initialise_bond_connectivity(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    // molecule chain, molecule k connects to k+1 for each component
    
    size_t num_molecules = p_parameters->num_part / p_parameters->chain_length;
    size_t num_bonds = (p_parameters->chain_length -1) * num_molecules;

    struct Bond *bonds =
        (struct Bond *)malloc(num_bonds * sizeof(struct Bond));

    /// \todo Specify bonds between particles, i.e., bonds[i].i and bonds[i].j for bonded particle pairs.

    size_t bond_index = 0;

    
    for (size_t molecule = 0; molecule < num_molecules; ++molecule)
    {
        size_t first = p_parameters->chain_length * molecule;

        for (size_t site = 0; site < p_parameters->chain_length -1; ++site)
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
    double r_eq = p_parameters->r_cut/(1+((2*p_parameters->kT)/(p_parameters->a[0]*p_parameters->r_cut)));
    size_t c_l = (size_t)p_parameters->chain_length;
    if (c_l == 0 || (double)c_l != p_parameters->chain_length
        || p_parameters->num_part % c_l != 0)
    {
        fprintf(stderr, "Particle count must be a multiple of a positive chain length.\n");
        exit(EXIT_FAILURE);
    }

    size_t num_molecules = p_parameters->num_part / c_l;
    size_t max_attempts = 100000;
    double min_separation = 0.5 * r_eq + 0.002;
    double min_separation_sq = min_separation * min_separation;
    double two_pi = 2.0 * acos(-1.0);
    struct Vec3D *candidate = malloc(c_l * sizeof(*candidate));
    if (candidate == NULL)
    {
        fprintf(stderr, "Could not allocate temporary chain positions.\n");
        exit(EXIT_FAILURE);
    }

    for (size_t molecule = 0; molecule < num_molecules; ++molecule)
    {
        size_t first = molecule * c_l;
        int type = p_vectors->type[first];
        int placed = 0;

        for (size_t attempt = 0; attempt < max_attempts && !placed; ++attempt)
        {
            double x_min = 0.0;
            double x_width = p_parameters->L.x;
            if (p_parameters->demix_initialisation)
            {
                x_width *= 0.5;
                if (type != 0)
                    x_min = x_width;
            }

            struct Vec3D position = v3(
                x_min + x_width * generate_uniform_random(),
                p_parameters->L.y * generate_uniform_random(),
                p_parameters->L.z * generate_uniform_random()
            );
            placed = 1;

            for (size_t site = 0; site < c_l; ++site)
            {
                struct Vec3D wrapped = v3_in_box(position, p_parameters->L);
                if (p_parameters->demix_initialisation
                    && ((type == 0 && wrapped.x >= 0.5 * p_parameters->L.x)
                        || (type != 0 && wrapped.x < 0.5 * p_parameters->L.x)))
                {
                    placed = 0;
                    break;
                }

                candidate[site] = wrapped;
                for (size_t previous = 0; previous < first; ++previous)
                {
                    struct Vec3D separation = v3_min_image(
                        v3_sub(wrapped, p_vectors->r[previous]),
                        p_parameters->L
                    );
                    if (v3_dot(separation, separation) < min_separation_sq)
                    {
                        placed = 0;
                        break;
                    }
                }
                if (!placed)
                    break;

                for (size_t previous_site = 0; previous_site < site; ++previous_site)
                {
                    struct Vec3D separation = v3_min_image(
                        v3_sub(wrapped, candidate[previous_site]),
                        p_parameters->L
                    );
                    if (v3_dot(separation, separation) < min_separation_sq)
                    {
                        placed = 0;
                        break;
                    }
                }
                if (!placed)
                    break;

                if (site + 1 < c_l)
                {
                    double cos_theta = 2.0 * generate_uniform_random() - 1.0;
                    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
                    double phi = two_pi * generate_uniform_random();
                    struct Vec3D step = v3(
                        r_eq * sin_theta * cos(phi),
                        r_eq * sin_theta * sin(phi),
                        r_eq * cos_theta
                    );
                    position = v3_add(position, step);
                }
            }
        }

        if (!placed)
        {
            fprintf(stderr,
                    "Could not place chain %zu without overlaps after %zu attempts.\n",
                    molecule,
                    max_attempts);
            free(candidate);
            exit(EXIT_FAILURE);
        }

        for (size_t site = 0; site < c_l; ++site)
            p_vectors->r[first + site] = candidate[site];
    }

    free(candidate);
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