#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "random.h"
#include "initialise.h"

// Initializes particle properties such as type, radius, and mass based on parameters.
void initialise_particles(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    int *type = p_vectors->type;
    double *radius = p_vectors->radius;
    double *mass = p_vectors->mass;
    double rho = (p_parameters->density);
    double R_max = p_parameters->R_max;
    double R_min = p_parameters->R_min;
    size_t num_part = p_parameters->num_part;
    for (size_t i = 0; i < num_part; i++)
    {
        type[i] = 0;
        radius[i] = R_min + (R_max-R_min)*generate_uniform_random();
        double V = PI*(4.0/3.0)*(radius[i]*radius[i]*radius[i]);
        mass[i] = V*rho;
    }
}

// Main initialization function that seeds the random generator, sets time to zero,
// and calls sub-functions to initialize particles, positions, and velocities.
void initialise(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    srand(SEED); // Positive integer as seed for random number generator
    p_vectors->time = 0.0; // Initialize the time to zero
    initialise_particles(p_parameters, p_vectors);
    if (p_parameters->run == RUN_COLLISION)
    {
        initialise_collision(p_parameters, p_vectors);
        return;
    }
    if (p_parameters->run == RUN_SILO)
        initialise_positions_silo(p_parameters, p_vectors);
    else
        initialise_positions(p_parameters, p_vectors);
    initialise_velocities(p_parameters, p_vectors);
    return;
}

// Part A: the test collision (Tasks A1 and A2).
void initialise_collision(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    /// \todo Task A1 and A2: set up the test collision. num_part is 1 or 2
    /// (NUM_TEST in setparameters.c) and gravity is off. Give the particles
    /// equal radii R_min and the matching mass (density * 4/3 pi R^3), and
    /// zero angular velocity. One particle (A1): place it just above the bottom
    /// wall, at height 1.05 * R in the middle of the box, moving straight
    /// down with speed v_impact. Two particles (A2): place them on a line
    /// along x in the middle of the box, centres 2 * 1.05 * R apart, moving
    /// towards each other with speed v_impact each. Until this is done the
    /// lattice of the starting code is used, at rest.
    // Equal radii R_min and the matching mass (this overrides the random radii
    // of initialise_particles), no rotation. The particles start 0.05 R apart
    // from touching (centre at 1.05 R from the wall, or 2.1 R from each other),
    // so the contact begins during the run and the tangential histories
    // (zeroed in alloc_nbrlist / alloc_colllist) start empty.
    const double R = p_parameters->R_min;
    const double m = p_parameters->density * (4.0 / 3.0) * PI * R * R * R;
    const double v0 = p_parameters->v_impact;
    const struct Vec3D mid = {0.5 * p_parameters->L.x, 0.5 * p_parameters->L.y, 0.5 * p_parameters->L.z};
    for (size_t i = 0; i < p_parameters->num_part; ++i)
    {
        p_vectors->radius[i] = R;
        p_vectors->mass[i] = m;
        p_vectors->v[i] = (struct Vec3D){0.0, 0.0, 0.0};
        p_vectors->omega[i] = (struct Vec3D){0.0, 0.0, 0.0};
    }
    if (p_parameters->num_part == 1) // A1: above the bottom wall, moving down
    {
        p_vectors->r[0] = (struct Vec3D){mid.x, mid.y, 1.05 * R};
        p_vectors->v[0].z = -v0;
    }
    else // A2: two particles on a line along x, moving towards each other
    {
        p_vectors->r[0] = (struct Vec3D){mid.x - 1.05 * R, mid.y, mid.z};
        p_vectors->r[1] = (struct Vec3D){mid.x + 1.05 * R, mid.y, mid.z};
        p_vectors->v[0].x = v0;
        p_vectors->v[1].x = -v0;
    }
}

// Initializes particle positions on a cubic lattice based on the box dimensions and particle count.
// Positions are evenly spaced, and the lattice size is adjusted to accommodate the largest particle.
// This is the bed of Part B, and the stand-in for the parts not yet written.
void initialise_positions(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    struct Vec3D dr;
    struct Index3D n;
    double dl;
    size_t ipart;
    size_t num_part = p_parameters->num_part;
    double *R = p_vectors->radius;
    double R_max = 0;
    for (ipart =0; ipart<num_part; ++ipart)
        R_max = (R[ipart]>R_max? R[ipart]: R_max);
    dl = 2.1*R_max; // lattice spacing just over one diameter of the largest particle, so no pair starts overlapping
    n.i = (int)floor(p_parameters->L.x / dl);
    n.j = (int)floor(p_parameters->L.y / dl);
    n.k = (int) (p_parameters->num_part/(n.i*n.j)+1);
    dr.x = p_parameters->L.x / (double)n.i;
    dr.y = p_parameters->L.y / (double)n.j;
    dr.z = dl;
    ipart = 0;
    for (size_t i = 0; i < n.i; ++i)
        for (size_t j = 0; j < n.j; ++j)
            for (size_t k = 0; k < n.k; ++k, ++ipart)
            {
                if (ipart >= num_part)
                    break;
                p_vectors->r[ipart].x = (i + 0.5) * dr.x;
                p_vectors->r[ipart].y = (j + 0.5) * dr.y;
                p_vectors->r[ipart].z = (k + 0.5) * dr.z;
            }
}

// Initializes the particle positions inside the silo (Task C3).
void initialise_positions_silo(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    /// \todo Task C3: fill the silo. Place the particles on a loose cubic
    /// lattice clipped to the cylinder, layer by layer from the bottom plate up:
    /// - spacing dl = 2.2 * R_max, so no two particles start in contact;
    /// - keep a site only when its distance to the axis (0.5*L.x, 0.5*L.y) is
    ///   below R_cyl - 1.1 * R_max, so no particle starts in the cylinder wall;
    /// - the first layer at height z_silo + 1.1 * R_max, just above the closed plate;
    /// - shift each coordinate of each site by a random amount of at most
    ///   0.04 * R_max (use generate_uniform_random), so the bottom layer does
    ///   not keep the square order of the lattice.
    /// For 7500 particles the lattice is about 100 mm high, against about 58 mm
    /// for the settled fill; the box has room for it. The fill then settles under gravity on the closed plate
    /// until t_open; no separate randomisation stage is needed. Until this is
    /// done the lattice of the starting code is used.
    // Square lattice of spacing dl, sites at axis + (i, j) dl, layers k from the
    // plate upwards; a site is kept when it lies inside the margin radius, and
    // each kept site gets a random shift of at most 0.04 R_max per coordinate.
    // The loops fill a whole layer before the next one starts.
    const double R_max = p_parameters->R_max;
    const double dl = 2.2 * R_max;
    const double R_site = p_parameters->R_cyl - 1.1 * R_max; // sites must lie closer to the axis than this
    const double jitter = 0.04 * R_max;
    const double x0 = 0.5 * p_parameters->L.x, y0 = 0.5 * p_parameters->L.y;
    const int m = (int)floor(R_site / dl);                    // sites at |i|, |j| <= m can be inside
    size_t ipart = 0;
    for (size_t k = 0; ipart < p_parameters->num_part; ++k)
    {
        double z = p_parameters->z_silo + 1.1 * R_max + k * dl;
        for (int i = -m; i <= m && ipart < p_parameters->num_part; ++i)
            for (int j = -m; j <= m && ipart < p_parameters->num_part; ++j)
            {
                double x = i * dl, y = j * dl;
                if (x * x + y * y >= R_site * R_site)
                    continue;
                p_vectors->r[ipart].x = x0 + x + jitter * (2.0 * generate_uniform_random() - 1.0);
                p_vectors->r[ipart].y = y0 + y + jitter * (2.0 * generate_uniform_random() - 1.0);
                p_vectors->r[ipart].z = z + jitter * (2.0 * generate_uniform_random() - 1.0);
                ++ipart;
            }
    }
}

void initialise_velocities(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    struct Vec3D sumv = (struct Vec3D){0.0, 0.0, 0.0};

    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        double sqrttgm = sqrt(p_parameters->Tg / p_vectors->mass[i]);
        p_vectors->v[i].x = sqrttgm * gauss();
        p_vectors->v[i].y = sqrttgm * gauss();
        p_vectors->v[i].z = sqrttgm * gauss();
        sumv.x += p_vectors->v[i].x;
        sumv.y += p_vectors->v[i].y;
        sumv.z += p_vectors->v[i].z;
    }

    sumv.x /= ((double)(p_parameters->num_part)); /* remove average velocity */
    sumv.y /= ((double)(p_parameters->num_part)); /* so total momentum is zero */
    sumv.z /= ((double)(p_parameters->num_part));
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->v[i].x -= sumv.x;
        p_vectors->v[i].y -= sumv.y;
        p_vectors->v[i].z -= sumv.z;
    }
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->omega[i] = (struct Vec3D){0.0, 0.0, 0.0};
    }
}
