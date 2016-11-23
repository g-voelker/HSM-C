/* computes the latitude dependent damping coefficient for the 
 * slab model after Pollard and Millard 1970 */
//#include <bits/mathcalls.h>
#include <math.h>
#include <stdio.h>
#include "../lib/constants.h"
#include "input.h"

double damping(double lat) {
  if (DBGFLG>2) { printf("  damping: setting local damping time scale\n"); fflush(NULL);}
  /* this is set to 1/5 days for North Atlantic
   * for world wide distribution see Park 2009 */
  double r0;
  /* set r0 in 1/s */
  r0  = 1.0/5.0/86400.0;
  return r0;

}

double coriolis(double lat) {
  if (DBGFLG>2) { printf("  coriolis: setting local Coriolis frequency\n"); fflush(NULL);}
  /* computes the Coriolis frequency at a given latitude */
  double f0, Om;
  Om  = 7.2921150e-5;
  f0  = 2.0*Om*sin(lat*2.0*PI/360.0);
  return f0;
}