// Independent checks of typing, mass updates and the neighbour-list result.
// Build with the normal sources except main.c; see code/README.md.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "structs.h"
#include "constants.h"
#include "setparameters.h"
#include "memory.h"
#include "initialise.h"
#include "fileoutput.h"
#include "nbrlist.h"
#include "forces.h"
#include "dynamics.h"

int main(int argc, char **argv)
{
    struct Parameters p;
    struct Vectors v;
    struct Nbrlist nl;
    set_parameters(&p);
    if (argc == 2)
        snprintf(p.restart_in_filename, sizeof(p.restart_in_filename), "%s", argv[1]);
    alloc_memory(&p, &v, &nl);
    load_restart(&p, &v);
    boundary_conditions(&p, &v);
    initialise_types(&p, &v);
    initialise_structure(&p, &v, &nl);
    build_nbrlist(&p, &v, &nl);
    if (v.num_bonds != 0) { fprintf(stderr, "This test is for B3's unbonded system.\n"); return 1; }

    const double E = calculate_forces(&p, &nl, &v);
    const size_t N = p.num_part;
    struct Vec3D *reference = calloc(N, sizeof(*reference));
    double Eref = 0.0, Wref = 0.0;
    size_t counts[2] = {0, 0}, pair_count = 0;
    size_t cutoff_crossings[3] = {0, 0, 0};
    const double scaling_steps[3] = {1e-6, 1e-7, 1e-8};
    double closest_cutoff_gap = 1e99;
    double M = 0.0;
    struct Vec3D momentum = {0, 0, 0};
    for (size_t i = 0; i < N; ++i)
    {
        const int expected_type = (i % 5 == 0 || i % 5 == 4) ? TYPE_CH3 : TYPE_CH2;
        if (v.type[i] != expected_type) return 2;
        counts[v.type[i]]++;
        M += p.mass[v.type[i]];
        momentum = v3_add(momentum, v3_scl(p.mass[v.type[i]], v.v[i]));
        for (size_t j = i + 1; j < N; ++j)
        {
            struct Vec3D d = v3_sub(v.r[i], v.r[j]);
            // Independent nearest-image calculation, valid for arbitrary wrapping.
            d.x -= p.L.x * nearbyint(d.x/p.L.x);
            d.y -= p.L.y * nearbyint(d.y/p.L.y);
            d.z -= p.L.z * nearbyint(d.z/p.L.z);
            const double r = v3_norm(d);
            closest_cutoff_gap = fmin(closest_cutoff_gap, fabs(r-p.r_cut));
            for (int q = 0; q < 3; ++q)
                if (r*(1-scaling_steps[q]) < p.r_cut && r*(1+scaling_steps[q]) > p.r_cut)
                    cutoff_crossings[q]++;
            if (r >= p.r_cut) continue;
            pair_count++;
            const double sig = (p.sigma[v.type[i]] + p.sigma[v.type[j]])/2;
            const double eps = sqrt(p.epsilon[v.type[i]] * p.epsilon[v.type[j]]);
            const double x = pow(sig/r, 6), xc = pow(sig/p.r_cut, 6);
            Eref += 4*eps*(x*x-x-xc*xc+xc);
            const double factor = 24*eps*(2*x*x-x)/(r*r);
            reference[i] = v3_add(reference[i], v3_scl(factor, d));
            reference[j] = v3_sub(reference[j], v3_scl(factor, d));
            Wref += factor*r*r;
        }
    }
    double max_force_error = 0;
    for (size_t i = 0; i < N; ++i)
        max_force_error = fmax(max_force_error, v3_norm(v3_sub(v.f[i], reference[i]))/fmax(1., v3_norm(reference[i])));
    const double energy_error = fabs(E-Eref)/fmax(1., fabs(Eref));
    const double W = 3*p.L.x*p.L.y*p.L.z*v.press_vir_nb;
    const double virial_error = fabs(W-Wref)/fmax(1., fabs(Wref));
    printf("types: CH3=%zu CH2=%zu; mass=%.12g amu; COM speed=%.12g internal\n", counts[0], counts[1], M, v3_norm(momentum)/M);
    printf("all-pairs reference: interacting pairs=%zu E=%.15g W=%.15g\n", pair_count, E, W);
    printf("all-pairs relative errors: energy=%.12g max force=%.12g virial=%.12g\n", energy_error, max_force_error, virial_error);
    printf("closest pair to cutoff: gap=%.12g angstrom\n", closest_cutoff_gap);
    printf("pairs crossing cutoff in virial scaling: h=1e-6: %zu; h=1e-7: %zu; h=1e-8: %zu\n",
           cutoff_crossings[0], cutoff_crossings[1], cutoff_crossings[2]);

    // A uniform test force must accelerate the two masses differently.
    const struct Vec3D force = {1., -2., 3.};
    for (size_t i = 0; i < N; ++i) { v.v[i] = (struct Vec3D){0,0,0}; v.f[i] = force; }
    const double K = update_velocities_half_dt(&p, &nl, &v);
    double Kref = 0, velocity_error = 0;
    for (size_t i = 0; i < N; ++i)
    {
        const double mass = (i % 5 == 0 || i % 5 == 4) ? 15.035 : 14.027;
        const struct Vec3D expected = v3_scl(p.dt/(2*mass), force);
        velocity_error = fmax(velocity_error, v3_norm(v3_sub(v.v[i], expected))/v3_norm(expected));
        Kref += 0.5*mass*v3_dot(expected, expected);
    }
    const double kinetic_error = fabs(K-Kref)/Kref;
    printf("mass-update relative errors: velocity=%.12g kinetic energy=%.12g\n", velocity_error, kinetic_error);
    free(reference);
    free_memory(&v, &nl);
    return !(energy_error < 1e-10 && max_force_error < 1e-10 && virial_error < 1e-10 && velocity_error < 1e-14 && kinetic_error < 1e-14);
}
