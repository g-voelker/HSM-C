/* computes the latitude dependent damping coefficient for the 
 * slab model after Pollard and Millard 1970 */
//#include <bits/mathcalls.h>
#include <math.h>
#include "../lib/constants.h"

double damping(double lat) {
  /* this is set to 1/5 days for North Atlantic
   * for world wide distribution see Park 2009 */
  double r0;
  /* set r0 in 1/s */
  r0  = 1/5/86400;
  return r0;
}

double coriolis(double lat) {
  /* computes the Coriolis frequency at a given latitude */
  double f0, Om;
  Om  = 7.2921150e-5;
  f0  = 2*Om*sin(lat*2*PI/360);
  return f0;
}