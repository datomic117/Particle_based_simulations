#include "structs.h"
#include "random.h"

// Generate a uniform random number between 0 and 1, inclusive
// Note: this is NOT the best random number out there, but it's quick and simple
double generate_uniform_random(void)
{
  double r;
  r = (double)rand() / (double)RAND_MAX;
  return r;
}



// Edited to get values between [-sqrt(3);sqrt(3)]
double gauss(void)
{
  return (generate_uniform_random() - 0.5) * (2.0 * sqrt(3.0));
}
