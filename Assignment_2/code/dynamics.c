#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// This function updates particle positions using their velocities.
// The positions are advanced by one full time step (dt), and displacement vectors
// (dr) for one time step are updated. The displacement since the last neighbor 
// list creation (stored in p_nbrlist->dr) is also updated for each particle.
void update_positions(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D dr_loc;
    struct Vec3D *r = p_vectors->r;     // Particle positions
    struct Vec3D *dr = p_vectors->dr;   // Displacement in one timestep
    struct Vec3D *v = p_vectors->v;     // Particle velocities
    struct DeltaR *dr_nbrlist = p_nbrlist->dr;  // Displacement since last neighbor list creation
    size_t num_part = p_parameters->num_part;
    double dt = p_parameters->dt;

    // Loop over all particles to update their positions
    for (size_t i = 0; i < num_part; i++)
    {
        dr_loc = v3_scl(dt, v[i]);   // Compute displacement in one timestep
        dr[i] = dr_loc;              // Store the displacement for this timestep
        r[i] = v3_add(r[i], dr_loc); // Update position

        // Track the accumulated displacement since the neighbor list was last
        // built; update_nbrlist uses it to decide when a rebuild is needed
        dr_nbrlist[i].v = v3_add(dr_nbrlist[i].v, dr_loc);
        dr_nbrlist[i].sq = v3_dot(dr_nbrlist[i].v, dr_nbrlist[i].v);
    }
}

// This function updates particle velocities by half a time step using the current forces.
// The updated velocities are used in the velocity-Verlet integration scheme.
// The function also calculates and returns the kinetic energy of the system.
double update_velocities_half_dt(struct Parameters *p_parameters,
                                 struct Nbrlist *p_nbrlist,
                                 struct Vectors *p_vectors)
{
    (void)p_nbrlist;

    double Ekin = 0.0;
    const double dt_half = 0.5 * p_parameters->dt;

    struct Vec3D *v = p_vectors->v;
    struct Vec3D *f = p_vectors->f;

    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        const int type = p_vectors->type[i];
        const double mass = p_parameters->mass[type];

        // Advance velocity by half a timestep using this site's mass.
        v[i] = v3_add(v[i], v3_scl(dt_half / mass, f[i]));

        // Accumulate kinetic energy using the same mass.
        Ekin += 0.5 * mass * v3_dot(v[i], v[i]);
    }

    const double volume =
        p_parameters->L.x *
        p_parameters->L.y *
        p_parameters->L.z;

    p_vectors->press_kin = 2.0 * Ekin / (3.0 * volume);

    return Ekin;
}
// This function applies periodic boundary conditions to ensure particles stay inside the simulation box.
// If a particle moves beyond the box, it is wrapped around to the opposite side.
// Keeping all positions inside [0,L) is what allows the cheap minimum image
// convention (v3_min_image) used everywhere else in the code.
void boundary_conditions(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    struct Vec3D invL;  // Inverse of the box size
    struct Vec3D *r = p_vectors->r;  // Particle positions
    struct Vec3D L = p_parameters->L;  // Box dimensions
    size_t num_part = p_parameters->num_part;  // Number of particles

    invL.x = 1.0 / L.x;
    invL.y = 1.0 / L.y;
    invL.z = 1.0 / L.z;

    // Loop over all particles and apply periodic boundary conditions
    for (size_t i = 0; i < num_part; i++)
    {
        r[i].x -= L.x * floor(r[i].x * invL.x);  // Apply periodic boundary in x-direction
        r[i].y -= L.y * floor(r[i].y * invL.y);  // Apply periodic boundary in y-direction
        r[i].z -= L.z * floor(r[i].z * invL.z);  // Apply periodic boundary in z-direction
    }
}

// This function applies a thermostat to maintain the system's temperature.
void thermostat(struct Parameters *p_parameters,
                struct Vectors *p_vectors,
                double Ekin)
{
    const double dof =
        3.0 * (double)p_parameters->num_part - 3.0;

    const double temperature = 2.0 * Ekin / dof;

    if (!isfinite(temperature) || temperature <= 0.0 ||
        p_parameters->tau_T <= 0.0)
    {
        fprintf(stderr,
                "Error: invalid thermostat temperature or coupling time.\n");
        exit(EXIT_FAILURE);
    }

    const double lambda_sq =
        1.0 + (p_parameters->dt / p_parameters->tau_T)
        * (p_parameters->kT / temperature - 1.0);

    if (!isfinite(lambda_sq) || lambda_sq <= 0.0)
    {
        fprintf(stderr, "Error: invalid thermostat scaling factor.\n");
        exit(EXIT_FAILURE);
    }

    const double lambda = sqrt(lambda_sq);

    for (size_t i = 0; i < p_parameters->num_part; i++)
    {
        p_vectors->v[i] = v3_scl(lambda, p_vectors->v[i]);
    }

    const double volume =
        p_parameters->L.x *
        p_parameters->L.y *
        p_parameters->L.z;

    p_vectors->press_kin =
        2.0 * Ekin * lambda_sq / (3.0 * volume);
}