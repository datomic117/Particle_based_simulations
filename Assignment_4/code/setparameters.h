#ifndef SETPARAMETERS_H_
#define SETPARAMETERS_H_
/* function  parameters */
/**
 * @brief Set the parameters struct used in the simulation. The set-up is
 * chosen with the macro RUN (see the top of setparameters.c), which can also be
 * given on the command line, e.g. gcc -O3 -DRUN=3 *.c -o dem -lm.
 *
 * @param[out] p_parameters Parameters of the simulation.
 */
void set_parameters(struct Parameters * p_parameters);

/**
 * @brief Sizes the silo and the box around it from the number of particles.
 *
 * The radius of the silo follows from the number of particles and the aspect
 * ratio H/R_cyl = 3 of the settled fill. The silo bottom plate sits at height
 * z_silo, one fill height above the floor, so that the discharged material has
 * room to fall and heap up below the orifice. The orifice diameter is D_ORIFICE
 * mean particle diameters. For 7500 particles this gives R_cyl = 17.7 mm
 * (a silo diameter of 17.7 d_p), z_silo = 53.1 mm and a box of
 * 70.8 x 70.8 x 185.9 mm.
 *
 * @param[in,out] p_parameters members used: num_part, R_min, R_max; set: R_cyl,
 * z_silo, R_orifice, orifice_open, L
 */
void set_silo_geometry(struct Parameters *p_parameters);

/**
 * @brief Applies the time-dependent part of a run, once per time step, before
 * the forces are calculated: in Part B it rotates gravity, in Part D it opens
 * the orifice.
 *
 * @param[in,out] p_parameters members used: run, t_tilt, tilt_rate, t_open; set: slope_angle, g, orifice_open
 * @param[in] time the current simulated time
 */
void update_parameters(struct Parameters *p_parameters, double time);

#endif /* SETPARAMETERS_H_ */
