#ifndef CONSTANTS_H_
#define CONSTANTS_H_

/* Define constants for the entire program here */

#define PI 3.141592653589
#define NUM_WALLS_MAX 10

/// Seed used for reproducible random initialization (can be changed for variability)
#define SEED 12345u
/// Magnitude of the gravitational acceleration (m/s^2)
#define GRAVITY 9.81

/// @name The three set-ups of assignment A4, selected by Parameters::run
///@{
#define RUN_COLLISION 1 //!< Part A: one particle hitting the bottom wall, or two particles colliding head-on
#define RUN_SLOPE 2     //!< Part B: a bed on the bottom wall of the box, tilted by rotating gravity
#define RUN_SILO 3      //!< Parts C and D: the silo, filled, settled and discharged through its orifice
///@}

#endif /* CONSTANTS_H_ */
