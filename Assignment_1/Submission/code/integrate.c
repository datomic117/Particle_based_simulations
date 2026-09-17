#include "integrate.h"
#include "forces.h"

/* File for the definition of the integrators and the integrator
switch from B5. */

// Advance one Euler-forward step using the current acceleration
double euler_step(
    int N,
    const double *m,
    Vec3D *r,
    Vec3D *v,
    Vec3D *a,
    double dt
) {
    // Update positions and velocities using the acceleration at t
    for (int i = 0; i < N; i++) {
        r[i] = v_add(r[i], V_scl(dt, v[i]));
        v[i] = v_add(v[i], V_scl(dt, a[i]));
    }

    // Recompute accelerations (and potential energy) at the new positions
    return compute_acc_and_potential(N, m, r, a);
}


// Advance one velocity-Verlet step using the current acceleration
double verlet_step(
    int N,
    const double *m,
    Vec3D *r,
    Vec3D *v,
    Vec3D *a,
    double dt
) {
    // Keep the accelerations at t for the velocity half-step below
    Vec3D *a_old = malloc(N * sizeof(Vec3D));

    for (int i = 0; i < N; i++) {
        a_old[i] = a[i];
    }

    // Update positions using the acceleration at t
    for (int i = 0; i < N; i++) {
        r[i] = v_add(
            r[i],
            v_add(V_scl(dt, v[i]), V_scl(0.5 * dt * dt, a_old[i]))
        );
    }

    // Recompute accelerations (and potential energy) at the new positions
    double U = compute_acc_and_potential(N, m, r, a);

    // Update velocities using the average of the accelerations at t and t+dt
    for (int i = 0; i < N; i++) {
        v[i] = v_add(
            v[i],
            V_scl(0.5 * dt, v_add(a_old[i], a[i]))
        );
    }

    free(a_old);

    return U;
}


// Runtime switch between the two integrators (B5)
double integrate_step(
    Integrator_type integrator,
    int N,
    const double *m,
    Vec3D *r,
    Vec3D *v,
    Vec3D *a,
    double dt
) {
    switch (integrator) {
        case EULER:
            return euler_step(N, m, r, v, a, dt);
        case VERLET:
            return verlet_step(N, m, r, v, a, dt);
        default:
            printf("Unknown integrator type\n");
            exit(1);
    }
}
