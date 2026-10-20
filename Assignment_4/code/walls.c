#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include "constants.h"
#include "structs.h"
#include "nbrlist.h"
#include "walls.h"

/**
 * @file walls.c
 * Wall interaction functions.
 *
 * A wall is a function with the signature
 *   bool wall(struct Parameters *p_parameters, double radius, struct Vec3D *r,
 *             struct DeltaR *riw, struct Vec3D *vw)
 * that returns true if the particle at centre position r with the given radius
 * overlaps the wall. When it returns true it must set
 *   riw -> the vector from the closest point on the wall to the particle centre,
 *          together with its squared length,
 *   vw  -> the velocity of the wall at that point (zero for a static wall).
 * A wall is activated by assigning it to p_parameters->wall_function[n] in
 * setparameters.c and counting it in p_parameters->num_walls.
 *
 * The six planar walls below bound the box. Which of them are active is decided
 * in setparameters.c; the directions without a wall are periodic.
 */

bool bottom_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw)
{
    (void)p_parameters;
    if (r->z < radius) //check if there is interaction with the bottom wall
    {
        double d = r->z;
        *riw = (struct DeltaR){0, 0, d, d * d};
        *vw = (struct Vec3D){0};
        return true;
    }
    else
    {
        return false;
    }
}

bool top_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw)
{
    if (r->z > p_parameters->L.z - radius) //check if there is interaction with the top wall
    {
        double d = r->z - p_parameters->L.z;
        *riw = (struct DeltaR){0, 0, d, d * d};
        *vw = (struct Vec3D){0};
        return true;
    }
    else
    {
        return false;
    }
}

bool left_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw)
{
    (void)p_parameters;
    if (r->x < radius) //check if there is interaction with the wall at x = 0
    {
        double d = r->x;
        *riw = (struct DeltaR){d, 0, 0, d * d};
        *vw = (struct Vec3D){0};
        return true;
    }
    else
    {
        return false;
    }
}

bool right_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw)
{
    if (r->x > p_parameters->L.x - radius) //check if there is interaction with the wall at x = Lx
    {
        double d = r->x - p_parameters->L.x;
        *riw = (struct DeltaR){d, 0, 0, d * d};
        *vw = (struct Vec3D){0};
        return true;
    }
    else
    {
        return false;
    }
}

bool front_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw)
{
    (void)p_parameters;
    if (r->y < radius) //check if there is interaction with the wall at y = 0
    {
        double d = r->y;
        *riw = (struct DeltaR){0, d, 0, d * d};
        *vw = (struct Vec3D){0};
        return true;
    }
    else
    {
        return false;
    }
}

bool back_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw)
{
    if (r->y > p_parameters->L.y - radius) //check if there is interaction with the wall at y = Ly
    {
        double d = r->y - p_parameters->L.y;
        *riw = (struct DeltaR){0, d, 0, d * d};
        *vw = (struct Vec3D){0};
        return true;
    }
    else
    {
        return false;
    }
}

/*
 * The silo. Both walls are centred on the vertical axis of the box, at
 * (0.5*L.x, 0.5*L.y). They are already registered in set_parameters for the
 * silo run (RUN_SILO); until you fill them in they never report a contact, so
 * the code compiles and runs, but nothing holds the particles in.
 */

/**
 * @brief A static vertical cylinder of radius R_cyl, centred on the axis of the
 * box, that forms the silo wall above the bottom plate (Task C1).
 */
bool cylinder_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw)
{
    /// \todo Task C1: the cylindrical silo wall. Its axis is the vertical line
    /// through (0.5*L.x, 0.5*L.y). With rc the distance of the particle centre
    /// to the axis, the particle overlaps the wall when
    /// rc > R_cyl - radius. The closest wall point is then at distance R_cyl
    /// from the axis, in the same direction as the particle, so riw is the
    /// horizontal vector (rc - R_cyl) * (dx, dy, 0) / rc, with (dx, dy) the
    /// horizontal offset of the centre from the axis, and riw.sq its square
    /// length. The wall is static: vw = 0. The cylinder only exists above the
    /// plate: return false for a particle with r->z < z_silo, so the discharged
    /// material can spread over the floor below. bottom_wall, at the top of
    /// walls.c, is an example of the contract.
    // Only above the plate; the cylinder is cut at z_silo.
    if (r->z < p_parameters->z_silo)
        return false;
    double dx = r->x - 0.5 * p_parameters->L.x; // horizontal offset of the centre from the axis
    double dy = r->y - 0.5 * p_parameters->L.y;
    double rc = sqrt(dx * dx + dy * dy);        // distance of the centre to the axis
    if (rc <= p_parameters->R_cyl - radius)
        return false; // rc > R_cyl - radius > 0 below, so rc = 0 (the axis) never gets here: no division by zero
    // The closest wall point is at R_cyl from the axis in the direction of the
    // particle, so riw = r - P = (rc - R_cyl) e_r, horizontal and pointing inwards.
    double d = rc - p_parameters->R_cyl;
    *riw = (struct DeltaR){d * dx / rc, d * dy / rc, 0.0, d * d};
    *vw = (struct Vec3D){0};
    return true;
}

/**
 * @brief The silo bottom: a static horizontal plate at height z_silo with a
 * circular orifice of radius R_orifice, centred on the axis, that is closed
 * while the silo is filled and settled and opened at t_open (Task C2).
 */
bool silo_bottom_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw)
{
    // While the orifice is closed (orifice_open = 0) the plate is a full
    // plane, which is the given part below: it works like bottom_wall, at
    // height z_silo. When the orifice is open, a particle whose centre is over
    // the hole falls through, unless it touches the edge of the hole.
    double dz = r->z - p_parameters->z_silo; // keep the sign: under the plate riw points down
    if (fabs(dz) >= radius)
        return false; // any contact, flat or edge, needs the plane within one radius

    double dx = r->x - 0.5 * p_parameters->L.x; // horizontal offset from the axis
    double dy = r->y - 0.5 * p_parameters->L.y;
    double rc = sqrt(dx * dx + dy * dy);         // distance of the centre to the axis

    if (!p_parameters->orifice_open || rc >= p_parameters->R_orifice)
    {
        // Contact with the flat part of the plate (or the closed plate).
        *riw = (struct DeltaR){0.0, 0.0, dz, dz * dz};
        *vw = (struct Vec3D){0};
        return true;
    }

    /// \todo Task C2: the edge of the hole. Here the orifice is open and the
    /// particle centre is over the hole (rc < R_orifice), within one radius of
    /// the plane of the plate. Such a particle can still touch the plate, at
    /// the edge of the hole. Work out which point of the plate is then closest
    /// to the centre, and return true when the particle overlaps it, with riw
    /// the vector from that point to the centre (riw.sq its square length) and
    /// vw = 0, as for every wall. Mind the particle that sits exactly on the
    /// axis. Until this is done, a particle over the hole never touches the
    /// plate: the silo discharges, but through a hole whose edge the particles
    /// do not feel.
    // The plate is the half-plane z = z_silo, rc' >= R_orifice. For a centre at
    // distance rc < R_orifice from the axis the horizontal distance to any plate
    // point (rc' >= R_orifice) is at least R_orifice - rc (triangle inequality),
    // with equality only on the rim circle in the azimuth of the particle:
    // P = axis + R_orifice e_r, at height z_silo. So riw = r - P =
    // ((rc - R_orifice) e_r, dz). On the axis every rim point is equally close:
    // take e_r = x-hat, the distance does not depend on the choice.
    double ex = 1.0, ey = 0.0; // e_r; the arbitrary direction on the axis
    if (rc > 0.0)
    {
        ex = dx / rc;
        ey = dy / rc;
    }
    double d = rc - p_parameters->R_orifice; // negative: the centre is inside the rim
    double sq = d * d + dz * dz;
    if (sq >= radius * radius)
        return false; // the rim is farther than one radius: free fall through the hole
    *riw = (struct DeltaR){d * ex, d * ey, dz, sq};
    *vw = (struct Vec3D){0};
    return true;
}

size_t count_in_silo(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    /// \todo Task D2: count the particles still in the silo, those whose centre
    /// is above the plate (r.z > z_silo). main prints the count as the last
    /// column of its table, which then gives N(t), the number of particles in
    /// the silo against time.
    // one loop over the particles: the centre above the plate means still in the silo
    size_t n = 0;
    for (size_t i = 0; i < p_parameters->num_part; ++i)
        if (p_vectors->r[i].z > p_parameters->z_silo)
            ++n;
    return n;
}
