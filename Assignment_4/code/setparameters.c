#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"
#include "walls.h"
#include "setparameters.h"

/*
 * Choosing a run. Every setting below has a default, and the ones you change
 * from run to run can also be given on the command line, so that every run is
 * one reproducible compile command, for example
 *
 *   gcc -O3 -DRUN=1 -DNUM_TEST=2 -DV_IMPACT=0.3 *.c -o dem -lm    Part A, head-on pair
 *   gcc -O3 -DRUN=2 *.c -o dem -lm                                 Part B, tilted bed
 *   gcc -O3 -DRUN=3 -DD_ORIFICE=8 *.c -o dem -lm                   Part D, orifice of 8 d_p
 *
 * RUN 1, 2 and 3 are RUN_COLLISION, RUN_SLOPE and RUN_SILO of constants.h.
 */
#ifndef RUN
#define RUN RUN_SLOPE      // which set-up: 1 collision (Part A), 2 slope (Part B), 3 silo (Parts C, D)
#endif
#ifndef NUM_TEST
#define NUM_TEST 1         // Part A: 1 particle against the bottom wall (A1), or 2 particles head-on (A2)
#endif
#ifndef V_IMPACT
#define V_IMPACT 0.3       // Part A: impact speed in m/s
#endif
#ifndef DT_FRACTION
#define DT_FRACTION 0.05   // the time step as a fraction of the contact time: 0.05 resolves a contact in 20 steps
#endif
#ifndef TILT_RATE
#define TILT_RATE 20.0     // Part B: rate at which the slope angle grows from zero, in degrees per second
#endif
#ifndef D_ORIFICE
#define D_ORIFICE 7.0      // Parts C, D: orifice diameter in mean particle diameters
#endif
#ifndef T_DISCHARGE
#define T_DISCHARGE 2.0    // Part D: simulated time after the orifice opens, in s
#endif

void set_silo_geometry(struct Parameters *p_parameters)
{
  double aspect = 3.0;                           // H/R_cyl of the settled fill, H of about 1.5 silo diameters
  double pack_frac = 0.6;                        // solids volume fraction of a settled random packing
  double R_mean = 0.5 * (p_parameters->R_min + p_parameters->R_max);
  double V_part = p_parameters->num_part * (4.0/3.0) * PI * R_mean * R_mean * R_mean;
  p_parameters->R_cyl = pow(V_part / (pack_frac * PI * aspect), 1.0/3.0);
  double H_fill = aspect * p_parameters->R_cyl;
  p_parameters->z_silo = H_fill;                 // catch region below the plate, one fill height deep
  p_parameters->R_orifice = 0.5 * D_ORIFICE * 2.0 * R_mean;
  p_parameters->orifice_open = 0;                // closed while the silo is filled and settles
  // Above the plate the box holds the loose starting lattice: clipped to the
  // cylinder it stacks to about 1.9 H_fill, and more
  // for a small silo, so 2.5 fill heights leave headroom.
  p_parameters->L = (struct Vec3D){4.0*p_parameters->R_cyl, 4.0*p_parameters->R_cyl,
                                   p_parameters->z_silo + 2.5*H_fill};
}

void set_parameters(struct Parameters *p_parameters)
/* Set the parameters of this simulation */
{
  p_parameters->run = RUN;
  p_parameters->density = 2500;        // mass density of the particle material, kg/m^3

  double R_min = 0.975e-3;
  double R_max = 1.05*R_min;
  p_parameters->R_min = R_min;         // minimum particle radius
  p_parameters->R_max = R_max;         // maximum particle radius; the mean diameter is 2.0 mm
  double kn = 1000;                    // spring stiffness of normal contact force
  double e_n_pp = 0.92, e_t_pp = 0.33, muf = 0.50;    // restitution and friction coefficients for particle-particle contacts
  double e_n_pw = 0.85, e_t_pw = 0.33, muf_w = 0.45;  // the same for particle-wall contacts, equal for all walls

  // The box of the starting code: periodic in x and y, a floor and a ceiling.
  // Parts A and B use it; the silo sizes its own box in set_silo_geometry.
  p_parameters->L = (struct Vec3D){4e1*R_min, 4e1*R_min, 1e2*R_max};
  p_parameters->num_walls = 2;
  p_parameters->wall_function[0] = bottom_wall;
  p_parameters->wall_function[1] = top_wall;
  p_parameters->g = (struct Vec3D){0.0, 0.0, -GRAVITY};

  // Fields of the other set-ups, given a value so that nothing is left unset.
  p_parameters->v_impact = 0.0;
  p_parameters->t_tilt = 0.0;
  p_parameters->tilt_rate = 0.0;
  p_parameters->slope_angle = 0.0;
  p_parameters->R_cyl = 0.0;
  p_parameters->z_silo = 0.0;
  p_parameters->R_orifice = 0.0;
  p_parameters->orifice_open = 0;
  p_parameters->t_open = 0.0;

  double t_end;                        // simulated time of the run; num_dt_steps follows once dt is known
  switch (p_parameters->run)
  {
  case RUN_COLLISION:
    // Part A: one or two particles, no gravity, so that the kinetic energy
    // before and after the collision gives the restitution coefficient.
    p_parameters->num_part = NUM_TEST;
    if (NUM_TEST != 1 && NUM_TEST != 2)
    {
      fprintf(stderr, "set_parameters: NUM_TEST must be 1 (wall) or 2 (pair), not %d\n", NUM_TEST);
      exit(EXIT_FAILURE);
    }
    p_parameters->v_impact = V_IMPACT;
    p_parameters->g = (struct Vec3D){0.0, 0.0, 0.0};
    t_end = 2e-3;                      // about ten contact times
    p_parameters->num_dt_printf = 1;   // every step, to resolve the collision
    p_parameters->num_dt_traj = 1000000000; // no trajectory needed
    break;
  case RUN_SLOPE:
    // Part B: a bed several layers thick settles under vertical gravity until
    // t_tilt, after which gravity rotates at a constant rate (update_parameters).
    p_parameters->num_part = 2000;
    p_parameters->t_tilt = 0.3;
    p_parameters->tilt_rate = TILT_RATE * PI / 180.0;
    t_end = p_parameters->t_tilt + 24.0 / TILT_RATE; // ramp from 0 to 24 degrees
    p_parameters->num_dt_printf = 100;
    p_parameters->num_dt_traj = 2000;  // a frame every 0.022 s
    break;
  case RUN_SILO:
    // Parts C and D: the fill settles on the closed plate until t_open, then
    // the orifice opens (update_parameters) and the silo discharges.
    p_parameters->num_part = 7500;
    set_silo_geometry(p_parameters);
    p_parameters->num_walls = 4;
    p_parameters->wall_function[2] = cylinder_wall;     // the silo wall (Task C1)
    p_parameters->wall_function[3] = silo_bottom_wall;  // the silo bottom plate with the orifice (Task C2)
    p_parameters->t_open = 0.3;
    t_end = p_parameters->t_open + T_DISCHARGE;
    p_parameters->num_dt_printf = 200; // a line every 2.2 ms
    p_parameters->num_dt_traj = 2000;  // a frame every 0.022 s
    break;
  default:
    fprintf(stderr, "set_parameters: unknown run %d\n", p_parameters->run);
    exit(EXIT_FAILURE);
  }

  // The box is at rest; set A and freq to oscillate it vertically instead.
  p_parameters->A = 0.0;
  p_parameters->freq = 0.0;

  strcpy(p_parameters->filename_xyz, "trajectories");  //filename (without extension) for the xyz file
  p_parameters->load_restart = 0;                      //if equal 1 restart file is loaded
  strcpy(p_parameters->restart_in_filename, "restart.dat");  //filename for loaded restart file
  p_parameters->num_dt_restart = 1000000000;                 // number of time steps between saves; the end state is always saved
  strcpy(p_parameters->restart_out_filename, "restart.dat"); //filename for saved restart file

  // During a contact the normal overlap behaves as a damped harmonic
  // oscillator, which can be solved analytically: half an oscillation period is
  // the contact time, and the damping decides how much of the approach velocity
  // remains at rebound, i.e. the restitution coefficient. The expressions below
  // run that solution backwards: from the spring constant kn and the prescribed
  // restitution coefficients they compute the contact time and the dashpot and
  // tangential-spring constants that reproduce those coefficients, for a
  // contact of two particles of the reference mass. (In the force routines the
  // dashpot constants are rescaled to the actual masses of the touching pair.)
  // A wall acts as a particle of infinite mass, so the mass of a particle-wall
  // contact is the particle mass itself rather than the reduced mass m/2 of a
  // pair; that is why the wall constants below differ from the pair constants.
  double mass_ref = p_parameters->density * (4.0 / 3.0) * PI * R_min * R_min * R_min; //mass_ref of a particle (later coefficients are corrected for real particle mass)
  double tcontact = sqrt(0.5 * mass_ref * (PI * PI + pow(log(e_n_pp), 2)) / kn); //computed contact time
  p_parameters->k_n_pp = kn;                                                     //normal elastic spring constant for particle-particle interactions
  p_parameters->eta_n_pp = -2.0 * log(e_n_pp) * (0.5 * mass_ref) / tcontact;     //normal dashpot damping constant for particle-particle interactions
  p_parameters->k_t_pp = (mass_ref / 7.0) * (PI * PI + pow(log(e_t_pp), 2)) / (tcontact * tcontact); //tangential elastic spring constant for particle-particle interactions
  p_parameters->eta_t_pp = -2.0 * log(e_t_pp) * (mass_ref / 7.0) / tcontact;     //tangential dashpot damping coeff. for particle-particle interactions
  p_parameters->fric_pp = muf;                                                   //friction coefficient for particle-particle interactions
  for(unsigned int i=0; i< p_parameters->num_walls; ++i)
  {
    p_parameters->k_n_pw[i] = mass_ref * (PI * PI + pow(log(e_n_pw), 2)) / (tcontact * tcontact);               //normal elastic spring constant for particle-wall interactions
    p_parameters->eta_n_pw[i] = -2.0 * log(e_n_pw) * (mass_ref) / tcontact;                                     //normal dashpot damping coeff. for particle-wall interactions
    p_parameters->k_t_pw[i] = (2.0 * mass_ref / 7.0) * (PI * PI + pow(log(e_t_pw), 2)) / (tcontact * tcontact); //tangential elastic spring constant for particle-wall interactions
    p_parameters->eta_t_pw[i] = -2.0 * log(e_t_pw) * (2.0 * mass_ref / 7.0) / tcontact;                         //tangential dashpot damping coeff. for particle-wall interactions
    p_parameters->fric_pw[i] = muf_w;
  }
  p_parameters->e_n_pp_set = e_n_pp;   // record what the contact parameters were tuned to
  for (unsigned int i = 0; i < p_parameters->num_walls; ++i)
    p_parameters->e_n_pw_set[i] = e_n_pw;
  p_parameters->mass_ref = mass_ref; //mass_ref of a particle (later coefficients are corrected for real particle mass)

  double v_small = 1e-2 * R_max/tcontact;          //a velocity scale
  p_parameters->Tg = 0.5 * mass_ref * v_small * v_small; //here Tg denotes the granular temperature, an average kinetic energy used for initialization
  p_parameters->r_cut = 2.0 * R_max;               //cut-off distance for pair interactions: the largest distance at which two particles can touch
  p_parameters->dt = DT_FRACTION * tcontact;       //integration time step: a contact must be resolved by many steps
  p_parameters->r_shell = 0.2 * p_parameters->r_cut;             //shell thickness for neighbor list
  p_parameters->num_dt_steps = (size_t)(t_end / p_parameters->dt + 0.5);

  if (p_parameters->r_cut > p_parameters->L.x / 2.0)
      fprintf(stderr, "Warning! r_cut > Lx/2\n");
  if (p_parameters->r_cut > p_parameters->L.y / 2.0)
    fprintf(stderr, "Warning! r_cut > Ly/2\n");
  if (p_parameters->r_cut > p_parameters->L.z / 2.0)
    fprintf(stderr, "Warning! r_cut > Lz/2\n");
}

void update_parameters(struct Parameters *p_parameters, double time)
{
  if (p_parameters->run == RUN_SLOPE)
  {
    /// \todo Task B: tilt the bed by rotating gravity. Before t_tilt the slope
    /// angle is zero (the bed settles); from t_tilt on it grows as
    /// theta = tilt_rate * (time - t_tilt). Store theta in
    /// p_parameters->slope_angle and set
    /// p_parameters->g = GRAVITY * (sin(theta), 0, -cos(theta)).
    // theta is zero up to t_tilt (g stays vertical) and then grows linearly;
    // the walls are not moved, only the direction of g in the box frame.
    double theta = 0.0;
    if (time > p_parameters->t_tilt)
      theta = p_parameters->tilt_rate * (time - p_parameters->t_tilt);
    p_parameters->slope_angle = theta;
    p_parameters->g = (struct Vec3D){GRAVITY * sin(theta), 0.0, -GRAVITY * cos(theta)};
  }
  if (p_parameters->run == RUN_SILO)
  {
    /// \todo Task C2: open the orifice at t_open: set p_parameters->orifice_open
    /// to 1 once time >= t_open. Your silo_bottom_wall reads this flag.
    // closed while the fill settles, open from t_open on (the flag is never reset)
    if (time >= p_parameters->t_open)
      p_parameters->orifice_open = 1;
  }
}
