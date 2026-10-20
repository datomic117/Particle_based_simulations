/*
 * Standalone check of the silo walls of walls.c (Tasks C1 and C2). It links the real walls.c; it is
 * not part of the program (it lives in a subdirectory, so `gcc *.c` in code/ does not see it).
 *
 *   gcc -O2 -Wall -Wextra -Icode code/tests/wall_check.c code/walls.c -o runs/_check/wall_check.exe -lm
 *   runs/_check/wall_check.exe > data/C_wall_check.txt
 *
 * The reference is brute force and does not use the formulas of walls.c: the closest point on the wall is
 * searched on a finely sampled circle (rim of the hole, or the line of the cylinder at the height of the
 * particle) with a golden-section refinement of the best sample. For random particles it compares
 * r_iw, riw.sq and the overlap flag, and reports the largest deviation in units of the particle radius.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "constants.h"
#include "structs.h"
#include "walls.h"

#define PI_ 3.14159265358979323846
#define M_COARSE 720 /* 0.5 degree; the refinement below finds the minimum between the neighbours of the best sample */

static double urand(void) { return (double)rand() / (double)RAND_MAX; }

/* squared distance from the centre (dx, dy, z) (offsets from the axis) to the circle point at angle th, radius Rw, height zp */
static double dist2(double dx, double dy, double z, double Rw, double zp, double th)
{
    double px = Rw * cos(th), py = Rw * sin(th);
    return (dx - px) * (dx - px) + (dy - py) * (dy - py) + (z - zp) * (z - zp);
}

/* closest point on the circle of radius Rw at height zp: coarse sampling, then golden-section search around the best sample */
static void brute_closest(double dx, double dy, double z, double Rw, double zp, double *px, double *py)
{
    int kbest = 0;
    double dbest = 1e300;
    for (int k = 0; k < M_COARSE; ++k)
    {
        double d = dist2(dx, dy, z, Rw, zp, 2.0 * PI_ * k / M_COARSE);
        if (d < dbest) { dbest = d; kbest = k; }
    }
    double h = 2.0 * PI_ / M_COARSE;
    double a = 2.0 * PI_ * kbest / M_COARSE - h, b = 2.0 * PI_ * kbest / M_COARSE + h;
    const double gr = 0.6180339887498949;
    double c = b - gr * (b - a), d = a + gr * (b - a);
    for (int it = 0; it < 60; ++it)
    {
        if (dist2(dx, dy, z, Rw, zp, c) < dist2(dx, dy, z, Rw, zp, d)) b = d; else a = c;
        c = b - gr * (b - a);
        d = a + gr * (b - a);
    }
    double th = 0.5 * (a + b);
    *px = Rw * cos(th);
    *py = Rw * sin(th);
}

struct Stats
{
    long n, n_overlap, n_flag_mismatch, n_borderline;
    double max_riw, max_sq; /* in units of the radius (riw) and radius^2 (sq) */
    long n_nan;
    long n_resolved;        /* overlaps where the direction of riw was compared directly */
    long n_notclosest;      /* the point of the code is off the rim or farther than the brute-force point */
};

/* Rw: radius of the circle of wall points. check_direction: compare the vector riw with the brute-force one directly. This is
 * only possible where the search can resolve the azimuth: very near the axis the distance to the rim hardly depends on the
 * azimuth (a shift of 1e-6 rad changes it by less than rounding), so there only |riw| and riw.z are compared. In all cases the
 * point of the code, P = r - riw, must lie on the circle and no sampled circle point may be nearer (that is what "closest" means). */
static void compare(struct Stats *s, bool flag, struct DeltaR riw, double dx, double dy, double z, double zp,
                    double px, double py, double radius, double Rw, bool check_direction)
{
    double bx = dx - px, by = dy - py, bz = z - zp;
    double bsq = bx * bx + by * by + bz * bz;
    bool bflag = bsq < radius * radius;
    s->n++;
    if (bflag) s->n_overlap++;
    if (fabs(sqrt(bsq) - radius) < 1e-9 * radius) { s->n_borderline++; return; } /* too close to the contact boundary to decide */
    if (flag != bflag) { s->n_flag_mismatch++; return; }
    if (!flag) return;
    if (!isfinite(riw.x) || !isfinite(riw.y) || !isfinite(riw.z) || !isfinite(riw.sq)) { s->n_nan++; return; }
    double dsq = fabs(riw.sq - bsq) / (radius * radius);
    if (dsq > s->max_sq) s->max_sq = dsq;
    double dz_ = fabs(riw.z - bz) / radius, dxy;
    if (check_direction) dxy = sqrt((riw.x - bx) * (riw.x - bx) + (riw.y - by) * (riw.y - by)) / radius;
    else dxy = fabs(sqrt(riw.x * riw.x + riw.y * riw.y) - sqrt(bx * bx + by * by)) / radius; /* on the axis only the length is defined */
    double dm = dz_ > dxy ? dz_ : dxy;
    if (dm > s->max_riw) s->max_riw = dm;
    if (check_direction) s->n_resolved++;
    double Px = dx - riw.x, Py = dy - riw.y, Pz = z - riw.z;
    if (fabs(sqrt(Px * Px + Py * Py) - Rw) > 1e-12 * radius || fabs(Pz - zp) > 1e-12 * radius) s->n_notclosest++;
    if (sqrt(riw.x * riw.x + riw.y * riw.y + riw.z * riw.z) > sqrt(bsq) + 1e-12 * radius) s->n_notclosest++;
}

static void report(const char *name, struct Stats *s)
{
    printf("%-50s n=%6ld overlap=%6ld flag mismatches=%ld (borderline skipped %ld) NaN/vw errors=%ld not-closest=%ld | max|riw-brute|=%.1e R (%ld compared as vector) max|sq-brute|=%.1e R^2\n",
           name, s->n, s->n_overlap, s->n_flag_mismatch, s->n_borderline, s->n_nan, s->n_notclosest, s->max_riw, s->n_resolved, s->max_sq);
}

int main(void)
{
    srand(2026);
    struct Parameters p = {0};
    p.R_cyl = 0.0176999;
    p.z_silo = 0.0530997;
    p.L = (struct Vec3D){4.0 * p.R_cyl, 4.0 * p.R_cyl, 0.2};
    const double R_min = 0.975e-3, R_max = 1.05 * R_min;
    const double cx = 0.5 * p.L.x, cy = 0.5 * p.L.y;
    struct Vec3D vw_;
    struct DeltaR riw;
    long total_bad = 0;
    printf("wall_check: brute-force comparison of walls.c (R_cyl %.7g m, z_silo %.7g m, radius %.4g..%.4g m)\n\n", p.R_cyl, p.z_silo, R_min, R_max);

    /* --- C2, edge of the hole, orifice open ---------------------------------------------------------------- */
    /* the D = 7 d_p orifice, and two holes smaller than a particle so that a centre near the axis touches the rim */
    const double R_orif[3] = {0.00699562, 0.9e-3, 0.5e-3};
    for (int io = 0; io < 3; ++io)
    {
        p.R_orifice = R_orif[io];
        p.orifice_open = 1;
        for (int mode = 0; mode < 2; ++mode) /* 0: random points near the rim, 1: random points near the axis */
        {
            struct Stats s = {0};
            const long N = 100000;
            for (long i = 0; i < N; ++i)
            {
                double radius = R_min + (R_max - R_min) * urand();
                double rc = mode == 0 ? p.R_orifice - 1.2 * radius * urand() : (p.R_orifice < 1.2 * radius ? p.R_orifice : 1.2 * radius) * urand();
                if (rc < 0.0) rc = 0.0;
                double th = 2.0 * PI_ * urand();
                double z = p.z_silo + (2.0 * urand() - 1.0) * 1.2 * radius; /* above and below the plane, beyond the contact range */
                struct Vec3D r = {cx + rc * cos(th), cy + rc * sin(th), z};
                bool flag = silo_bottom_wall(&p, radius, &r, &riw, &vw_);
                double ddx = r.x - cx, ddy = r.y - cy; /* offsets as the code sees them (rounding of the box centre) */
                double px, py;
                brute_closest(ddx, ddy, z - p.z_silo, p.R_orifice, 0.0, &px, &py);
                compare(&s, flag, riw, ddx, ddy, z - p.z_silo, 0.0, px, py, radius, p.R_orifice, rc >= 0.1 * p.R_orifice);
                if (flag && (vw_.x != 0.0 || vw_.y != 0.0 || vw_.z != 0.0)) s.n_nan++; /* the wall is static */
            }
            char name[100];
            snprintf(name, sizeof name, "C2 edge, R_orifice=%.4g mm, %s", 1e3 * p.R_orifice, mode == 0 ? "near the rim" : "near the axis");
            report(name, &s);
            total_bad += s.n_flag_mismatch + s.n_nan + s.n_notclosest + (s.max_riw > 1e-6) + (s.max_sq > 1e-6);
        }
    }

    /* --- C2, on the axis exactly (rc = 0) ------------------------------------------------------------------ */
    {
        struct Stats s = {0};
        p.R_orifice = 0.5e-3;
        p.orifice_open = 1;
        for (int i = 0; i < 2000; ++i)
        {
            double radius = R_min + (R_max - R_min) * urand();
            double z = p.z_silo + (2.0 * urand() - 1.0) * 1.2 * radius;
            struct Vec3D r = {cx, cy, z};
            bool flag = silo_bottom_wall(&p, radius, &r, &riw, &vw_);
            double px, py;
            brute_closest(r.x - cx, r.y - cy, z - p.z_silo, p.R_orifice, 0.0, &px, &py);
            compare(&s, flag, riw, r.x - cx, r.y - cy, z - p.z_silo, 0.0, px, py, radius, p.R_orifice, false);
        }
        report("C2 exactly on the axis (rc = 0), R_orifice=0.5 mm", &s);
        total_bad += s.n_flag_mismatch + s.n_nan + s.n_notclosest + (s.max_riw > 1e-6) + (s.max_sq > 1e-6) + (s.n_overlap == 0);
        p.R_orifice = 0.00699562;
        struct Vec3D r = {cx, cy, p.z_silo + 0.3e-3};
        bool flag = silo_bottom_wall(&p, R_min, &r, &riw, &vw_);
        printf("   real orifice (R_orifice=%.4g mm), centre on the axis 0.3 mm above the plane: overlap=%d (the rim is %.2f mm away)\n",
               1e3 * p.R_orifice, (int)flag, 1e3 * p.R_orifice);
        if (flag) total_bad++;
    }

    /* --- C2, continuity at rc = R_orifice, and closed/open orientation ------------------------------------------ */
    printf("\ncontinuity at the rim (R_orifice = 6.996 mm, radius 1 mm, dz = 0.4 mm): edge branch at rc = R_orifice - eps against the flat plate (0, 0, dz) at rc = R_orifice\n");
    p.R_orifice = 0.00699562;
    p.orifice_open = 1;
    for (int e = 1; e <= 9; e += 2)
    {
        double eps = p.R_orifice * pow(10.0, -e);
        double radius = 1.0e-3, dz = 0.4e-3;
        struct Vec3D r = {cx + p.R_orifice - eps, cy, p.z_silo + dz};
        struct DeltaR edge, flat;
        bool f1 = silo_bottom_wall(&p, radius, &r, &edge, &vw_);
        r.x = cx + p.R_orifice; /* on the rim circle: flat branch (rc >= R_orifice) */
        bool f2 = silo_bottom_wall(&p, radius, &r, &flat, &vw_);
        double jump = sqrt((edge.x - flat.x) * (edge.x - flat.x) + (edge.y - flat.y) * (edge.y - flat.y) + (edge.z - flat.z) * (edge.z - flat.z));
        printf("   eps = 1e-%d R_orifice = %.3e m: flags %d %d, |riw_edge - riw_flat| = %.3e m, |sq_edge - sq_flat| = %.2e m^2\n",
               e, eps, (int)f1, (int)f2, jump, fabs(edge.sq - flat.sq));
        /* the jump must be eps (the shift of the centre) up to rounding of the position, ~1e-16 m */
        if (!f1 || !f2 || fabs(jump - eps) > 1e-12) total_bad++;
    }
    {
        long n_over = 0, bad_closed = 0, bad_open_far = 0, n_far = 0;
        for (long i = 0; i < 100000; ++i)
        {
            double radius = R_min + (R_max - R_min) * urand();
            double rc = 0.02 * urand(), th = 2.0 * PI_ * urand();
            double z = p.z_silo + (2.0 * urand() - 1.0) * 1.2 * radius;
            struct Vec3D r = {cx + rc * cos(th), cy + rc * sin(th), z};
            double dz = z - p.z_silo;
            bool expect = fabs(dz) < radius;
            if (expect) n_over++;
            p.orifice_open = 0; /* closed plate: a full plane everywhere, also over the hole */
            bool f = silo_bottom_wall(&p, radius, &r, &riw, &vw_);
            if (f != expect || (f && (riw.x != 0.0 || riw.y != 0.0 || fabs(riw.z - dz) > 1e-15 || fabs(riw.sq - dz * dz) > 1e-18))) bad_closed++;
            p.orifice_open = 1; /* open, outside the hole (rc >= R_orifice): the same flat plate */
            if (sqrt((r.x - cx) * (r.x - cx) + (r.y - cy) * (r.y - cy)) >= p.R_orifice)
            {
                f = silo_bottom_wall(&p, radius, &r, &riw, &vw_);
                if (f != expect || (f && (riw.x != 0.0 || riw.y != 0.0 || fabs(riw.z - dz) > 1e-15))) bad_open_far++;
                n_far++;
            }
        }
        printf("closed plate, 1e5 random points (rc up to 20 mm, also over the hole), %ld in contact: deviations from the flat plane = %ld\n", n_over, bad_closed);
        printf("open plate, %ld random points with rc >= R_orifice: deviations from the flat plane = %ld\n", n_far, bad_open_far);
        total_bad += bad_closed + bad_open_far;
    }

    /* --- C1, cylinder -------------------------------------------------------------------------------------------- */
    printf("\n");
    {
        struct Stats s = {0};
        long below = 0, bad_below = 0, inside_bad = 0;
        for (long i = 0; i < 100000; ++i)
        {
            double radius = R_min + (R_max - R_min) * urand();
            double rc = p.R_cyl - 1.2 * radius + 1.4 * radius * urand(); /* from outside the contact range to 0.2 R beyond the wall */
            double th = 2.0 * PI_ * urand();
            double z = p.z_silo + (urand() * 1.2 - 0.2) * 0.02; /* mostly above the plate, 1/6 below */
            struct Vec3D r = {cx + rc * cos(th), cy + rc * sin(th), z};
            bool flag = cylinder_wall(&p, radius, &r, &riw, &vw_);
            if (z < p.z_silo)
            {
                below++;
                if (flag) bad_below++;
                continue;
            }
            double ddx = r.x - cx, ddy = r.y - cy, px, py;
            brute_closest(ddx, ddy, z, p.R_cyl, z, &px, &py); /* the cylinder has a wall point at every height: zp = z */
            compare(&s, flag, riw, ddx, ddy, z, z, px, py, radius, p.R_cyl, true);
            if (flag && (vw_.x != 0.0 || vw_.y != 0.0 || vw_.z != 0.0 || riw.z != 0.0)) inside_bad++;
        }
        report("C1 cylinder, random points near the wall", &s);
        printf("   below the plate (z < z_silo): %ld points, %ld reported as contact; static wall (vw = 0, riw.z = 0) violated: %ld\n", below, bad_below, inside_bad);
        total_bad += s.n_flag_mismatch + s.n_nan + s.n_notclosest + (s.max_riw > 1e-6) + (s.max_sq > 1e-6) + bad_below + inside_bad;
        struct Vec3D r = {cx, cy, p.z_silo + 0.01};
        bool flag = cylinder_wall(&p, R_max, &r, &riw, &vw_);
        printf("   particle exactly on the axis: overlap=%d (rc = 0 < R_cyl - radius, so the division by rc is never reached)\n", (int)flag);
        if (flag) total_bad++;
    }

    printf("\n%s (%ld failed checks)\n", total_bad == 0 ? "ALL CHECKS PASSED" : "CHECKS FAILED", total_bad);
    return total_bad == 0 ? 0 : 1;
}
