#include <stdio.h>
#include <math.h>
//#include <netcdf.h>
#include "complex.h"
#include "alloc_space.h"
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



int main(void) {
  if (DBGFLG>1) {printf("main: initialize variables\n"); fflush(NULL);}
  // iterator integers
  int nn, mm;
  
  // land sea mask related variables
  char lsmfile[] = "lsm-hres.nc";
  int lsmid, latid, lonid, ncid
  size_t NLAT, NLON;
  double ***mask, *lat, *lon;
  
  // subsetting related variables
  int NLATMIN, NLATMAX, NLONMIN, NLONMAX;
  double minimum, maximum;

  // time loop related variables
  double *taux, *tauy, *time, *mld;
  int leap;

  if (DBGFLG>1) {printf("main: set land sea mask\n"); fflush(NULL);}
  lsm(lsmfile, ncid, lsmid, latid, lonid, NLAT, NLON, mask, lat, lon);
  
  if (DBGFLG>1) {printf("main: get subset\n"); fflush(NULL);}
  for (nn=0; nn<NLAT; nn++){ // iterate over all latitudes
    if (minimum>abs(lat[nn]-LATMIN){
      NLATMIN = nn;
      minimum = lat[nn];
    }
    if (maximum>abs(lat[nn]-LATMAX){
      NLATMAX = nn;
      maximum = lat[nn];
    }
  }
  for (nn=0; nn<NLON; nn++){ // iterate over all longitudes
    if (minimum>abs(lon[nn]-LONMIN){
      NLONMIN = nn;
      minimum = lon[nn];
    }
    if (maximum>abs(lon[nn]-LONMAX){
      NLONMAX = nn;
      maximum = lon[nn];
    }
  }
  
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
  time  = dvector(0,8760 + leap*24);
  mld   = dvector(0,8760 + leap*24);
  taux  = dvector(0,8760 + leap*24);
  tauy  = dvector(0,8760 + leap*24);
  
  // note: latitudes are sorted the wrong way around.
  for (nn=NLATMAX; nn<=NLATMIN; nn++){
    for (mm=NLONMIN; mm<=NLONMAX; mm++){
      getdata(nn, mm, time, mld, taux, tauy)
      // get mld
      // interpolate mld
      // calculate u and v
      // write out mld, time, u and v
    }
  }









}
