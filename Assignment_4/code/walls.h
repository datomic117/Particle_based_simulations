#ifndef WALLS_H_
#define WALLS_H_

#include <stddef.h>
#include <stdbool.h>

/**
 * @brief Model for the wall located at z=0, the floor of the box that catches the discharged material.
 *
 * @param[in] p_parameters unused for this static wall
 * @param[in] radius of the particle
 * @param[in] r position vector of the particle
 * @param[out] riw position of particle minus closest point on wall rw: riw = r-rw
 * @param[out] vw velocity of the wall at point rw, zero for this static wall
 * @return bool, true if particle and wall overlap, false otherwise
 */
bool bottom_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw);

/**
 * @brief Model for the wall located at z=L.z.
 *
 * @param[in] p_parameters member used: L
 * @param[in] radius of the particle
 * @param[in] r position vector of the particle
 * @param[out] riw position of particle minus closest point on wall rw: riw = r-rw
 * @param[out] vw velocity of the wall at point rw, zero for this static wall
 * @return bool, true if particle and wall overlap, false otherwise
 */
bool top_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw);

/**
 * @brief Models for the four vertical side walls, at x=0, x=L.x, y=0 and y=L.y.
 * Activate them in set_parameters to replace the periodic boundaries in those
 * directions, for instance to confine a bed in a rectangular box.
 *
 * @param[in] p_parameters member used: L
 * @param[in] radius of a particle
 * @param[in] r position vector of a particle
 * @param[out] riw position of particle minus closest point on wall rw: riw = r-rw
 * @param[out] vw velocity of wall at point rw, zero for these static walls
 * @return bool, true if particle and wall overlap false otherwise
 */
bool left_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw);
bool right_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw);
bool front_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw);
bool back_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw);

/**
 * @brief Model for a static vertical cylinder of radius R_cyl, centred on the
 * axis of the box, that forms the silo wall above the bottom plate (Task C1).
 *
 * @param[in] p_parameters members used: L, R_cyl, z_silo
 * @param[in] radius of a particle
 * @param[in] r position vector of a particle
 * @param[out] riw position of particle minus closest point on wall rw: riw = r-rw
 * @param[out] vw velocity of wall at point rw, zero for this static wall
 * @return bool, true if particle and wall overlap false otherwise
 */
bool cylinder_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw);
/**
 * @brief Model for the silo bottom: a static horizontal plate at height z_silo
 * with a circular orifice of radius R_orifice, centred on the axis, that is
 * closed while the silo is filled and settles and opened at t_open (Task C2).
 *
 * @param[in] p_parameters members used: L, z_silo, R_orifice, orifice_open
 * @param[in] radius of a particle
 * @param[in] r position vector of a particle
 * @param[out] riw position of particle minus closest point on wall rw: riw = r-rw
 * @param[out] vw velocity of wall at point rw, zero for this static wall
 * @return bool, true if particle and wall overlap false otherwise
 */
bool silo_bottom_wall(struct Parameters *p_parameters, double radius, struct Vec3D *r, struct DeltaR *riw, struct Vec3D *vw);
/**
 * @brief The number of particles still in the silo: those whose centre is above
 * the bottom plate (Task D2). Printed by main as the last column of its table.
 *
 * @param[in] p_parameters members used: num_part, z_silo
 * @param[in] p_vectors member used: r
 * @return the number of particles with r.z > z_silo
 */
size_t count_in_silo(struct Parameters *p_parameters, struct Vectors *p_vectors);
#endif  /* WALLS_H_ */
