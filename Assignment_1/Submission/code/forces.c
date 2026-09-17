#include "forces.h"
#define G 6.6743015e-11

double compute_acc_and_potential(
    int N,
    const double *m,
    const Vec3D *r,
    Vec3D *a
) {
    // Set all accelerations to zero
    for (int i = 0; i < N; i++) {
        a[i] = v_out(0.0, 0.0, 0.0);
    }

    // Total potential energy
    double U = 0.0;

    // Loop over every unique pair
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {

            // Vector from body i to body j
            Vec3D rij_vec = v_sub(r[j], r[i]);

            // Distance between the bodies
            double rij = v_norm(rij_vec);

            // Common factor G/r^3
            double factor = G / (rij * rij * rij);

            // Acceleration of body i due to body j
            a[i] = v_add(
                a[i],
                v_scl(factor * m[j], rij_vec)
            );

            // Acceleration of body j due to body i
            a[j] = v_add(
                a[j],
                v_scl(-factor * m[i], rij_vec)
            );

            // Add this pair's potential energy
            U += -G * m[i] * m[j] / rij;
        }
    }

    return U;
}