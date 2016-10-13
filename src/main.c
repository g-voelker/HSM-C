#include <stdio.h>
#include <math.h>
#include "../lib/complex.h"
#include "../lib/alloc_space.h"
//#include "constants.h"
//#include "macros.h"
//#include "damping.h"
#include "lsm.h"
#include "input.h"

/*
void slab(int year, int nlat, int nlon) {
  // netcdf ids, land-sea-mask file
  int ncid, latid, lonid, lsmid;
  char lsm[40];

  // input data arays, parameters, output variables
  double *tauX, *tauY, *mld, *time;
  double lat, lon, ff, rr, rho0;
  double *u, *v;

  // compute damping and coriolis coefficients
  ff  = corioli(lat);
  rr  = damping(lat);

  // open land sea mask
}*/

struct dat2D {
  int nlat, nlon;
  double *lat, *lon;
  double **data;
};

int main(void) {
  if (DBGFLG > 1) {
    printf("main: initialize variables\n");
    fflush(NULL);
  }
  // iterator integers
  int nn, mm;

  // land sea mask related variables
  char lsmfile[] = "static/lsm-hres.nc";

  // subsetting related variables
  int NLATMIN, NLATMAX, NLONMIN, NLONMAX;
  double minimum, maximum;

  // time loop related variables
  double *taux, *tauy, *time, *mld;
  int leap;

  if (DBGFLG > 1) {
    printf("main: set land sea mask\n");
    fflush(NULL);
  }
  struct dat2D lsmask = lsm(lsmfile);

  if (DBGFLG>1) {printf("main: get subset\n"); fflush(NULL);}
  NLATMIN = NLATMAX = NLONMIN = NLONMAX = 0;
  minimum = maximum = 180;
  for (nn=0; nn<lsmask.nlat; nn++){ // iterate over all latitudes
    if (minimum>fabs(lsmask.lat[nn]-LATMIN)) {
      NLATMIN = nn;
      minimum = fabs(lsmask.lat[nn]-LATMIN);
    }
    if (maximum>fabs(lsmask.lat[nn]-LATMAX)){
      NLATMAX = nn;
      maximum = fabs(lsmask.lat[nn]-LATMAX);
    }
  }
  minimum = maximum = 360;
  for (nn=0; nn<lsmask.nlon; nn++){ // iterate over all longitudes
    if (minimum>fabs(lsmask.lon[nn]-LONMIN)){
      NLONMIN = nn;
      minimum = fabs(lsmask.lon[nn]-LONMIN);
    }
    if (maximum>fabs(lsmask.lon[nn]-LONMAX)){
      NLONMAX = nn;
      maximum = fabs(lsmask.lon[nn]-LONMAX);
    }
  }

  printf("%.2d; %.2d; %.2d; %.2d;\n", NLATMIN, NLATMAX, NLONMIN, NLONMAX);
  printf("%.2f; %.2f; %.2f; %.2f;\n", lsmask.lat[NLATMIN], lsmask.lat[NLATMAX], lsmask.lon[NLONMIN], lsmask.lon[NLONMAX]);

  if (DBGFLG>1) {printf("main: begin loop over points\n"); fflush(NULL);}
  // set pointers for loop
  // check for leap years
  if (YEAR%4==0){
    if (YEAR%100!=0){
      leap  = 1;
    } else {
      leap  = 0;
    }
  } else {
    leap  = 0;
  }
  // allocate pointers
  time  = dvector(0, 8760 + leap * 24);
  mld   = dvector(0, 8760 + leap * 24);
  taux  = dvector(0, 8760 + leap * 24);
  tauy  = dvector(0, 8760 + leap * 24);

  fflush(NULL);

  for (nn=NLATMIN; nn<=NLATMAX; nn++){
    for (mm=NLONMIN; mm<=NLONMAX; mm++){

      // getdata(nn, mm, time, mld, taux, tauy);
      // get mld
      // interpolate mld
      // calculate u and v
      // write out mld, time, u and v
    }
  }


  return -1;
}
