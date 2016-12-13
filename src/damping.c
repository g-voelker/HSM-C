/* computes the latitude dependent damping coefficient for the 
 * slab model after Pollard and Millard 1970 */
//#include <bits/mathcalls.h>
#include <math.h>
#include <stdio.h>
#include "../lib/constants.h"
#include "../lib/structs.h"
#include "input.h"

double damping(dat2d_2 *damp, int nlat, double lon) {
  if (DBGFLG>2) { printf("  damping: setting local damping time scale\n"); fflush(NULL);}
  double r0, aux;

  // use negative representation for NA
  if (lon>=180){
    aux = lon-360;
  } else {
    aux = lon;
  }
  // check if in NA and set r0 accordingly
  if (damp->xx[nlat]<=30.0){
    if ((aux>fmax((-27 / 15 * damp->xx[nlat] - 64.0), -100.0))&(aux<=25.0)){
      r0 = damp->y2[nlat];
    } else {
      r0 = damp->y1[nlat];
    }
  } else {
    if ((aux>-100.0)&(aux<=25.0)){
      r0 = damp->y2[nlat];
    } else{
      r0 = damp->y1[nlat];
    }
  }

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