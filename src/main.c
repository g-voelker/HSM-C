#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
//#include "constants.h"
//#include "macros.h"
#include "../lib/constants.h"
#include "../lib/dalloc.h"
#include "damping.h"
#include "solveode.h"
#include "getdata.h"
#include "header.h"
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

struct dat1D {
  int ntime;
  double *time;
  double *data;
};

struct dat2D {
  int nlat, nlon;
  double *lat, *lon;
  double **data;
};

struct dat3D {
  int nt, nlat, nlon;
  double *lat, *lon, *time;
  double ***data;
};

int main(void) {
  if (DBGFLG > 1) {
    printf("main: initialize variables\n");
    fflush(NULL);
  }
  // iterator integers
  int nn, mm;

  // constants
  double rho0 = RHO, r0, f0;

  // land sea mask related variables
  char lsmfile[] = "../static/lsm-hres.nc";

  // subsetting related variables
  int NLATMIN, NLATMAX, NLONMIN, NLONMAX;
  double minimum, maximum;

  // time loop related variables
  int leap;
  struct tm tcon;

  // FFT related variables
  double *freqs;
  fftw_complex *aux, *AUX;
  fftw_plan fft, ifft;

  if (DBGFLG > 2) {printf("main: set time axis\n");fflush(NULL);}

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

  // allocate variables and set time
  double taux[8760 + leap * 24], tauy[8760 + leap * 24],
          time[8760 + leap * 24], mld[8760 + leap * 24];

  for (nn=0; nn< (8760 + leap * 24); nn++){
    mld[nn] = taux[nn] = tauy[nn] = time[nn] = 0.0;
  }

  // set basic time constructor
  tcon.tm_year = YEAR - 1900;
  tcon.tm_mon = 0;
  tcon.tm_mday = 0;
  tcon.tm_wday = 0;
  tcon.tm_yday = 0;
  tcon.tm_hour = 0;
  tcon.tm_min = 0;
  tcon.tm_sec = 0;
  tcon.tm_gmtoff = 0;
  tcon.tm_isdst = 0;

  for (nn=0; nn<365+leap; nn++){
    tcon.tm_yday = nn;
    for (mm=0; mm<24; mm++){
      tcon.tm_hour = mm;
      time[nn*24 + mm] = mktime(&tcon);
    }
  }

  if (DBGFLG > 1) {printf("main: set land sea mask\n");fflush(NULL);}
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

  // prepare FFTs
  // set frequencies
  aux = fftw_malloc(sizeof(fftw_complex) * DFFT_LEN);
  AUX = fftw_malloc(sizeof(fftw_complex) * DFFT_LEN);
  freqs = (double*) fftw_malloc(sizeof(double) * DFFT_LEN);
  for (nn=0; nn<DFFT_LEN/2; nn++) {
    freqs[nn] = nn / DFFT_LEN / 3600 * 2 * PI;
  }
  for (nn=DFFT_LEN/2; nn<DFFT_LEN; nn++) {
    freqs[nn] = (DFFT_LEN - nn) / DFFT_LEN / 3600 * 2 * PI;
  }
  if (DBGFLG>2) {printf("main: set fftw plans\n"); fflush(NULL);}
  // define transforms
  fft = fftw_plan_dft_1d(DFFT_LEN, aux, AUX, FFTW_FORWARD, FFTW_ESTIMATE);
  ifft = fftw_plan_dft_1d(DFFT_LEN, AUX, aux, FFTW_BACKWARD, FFTW_ESTIMATE);

  if (DBGFLG>1) {
    printf("main: subset boundaries are set by:\n");
    printf("     (latmin, latmax, lonmin, lonmax):\n");
    printf("      %.2f; %.2f; %.2f; %.2f;\n",
           lsmask.lat[NLATMIN], lsmask.lat[NLATMAX], lsmask.lon[NLONMIN], lsmask.lon[NLONMAX]);
    printf("main: begin loop over points\n");
    fflush(NULL);
  }

  fflush(NULL);

  for (nn=NLATMIN; nn<=NLATMAX; nn++){
//  for (nn=NLATMIN; nn<NLATMIN+1; nn++){
    for (mm=NLONMIN; mm<=NLONMAX; mm++){
//    for (mm=NLONMIN; mm<NLONMIN+1; mm++){
      if (DBGFLG>2) { printf("    (%d, %d)\n", nn, mm); fflush(NULL);}
      // load stress and mixed layer depth data
      getdata(nn, mm, leap, time, mld, taux, tauy);
      // set damping parameter
      r0 = damping(lsmask.lat[nn]);
      // set Coriolis frequency
      f0 = coriolis(lsmask.lat[nn]);
      // calculate u and v
      solve(fft, ifft, r0, f0, rho0, leap, taux, tauy, mld, freqs, aux, AUX);
      // write out mld, time, u and v
    }
  }

  if (DBGFLG>2) {printf("main: clean up slab model\n");fflush(NULL);}

  fftw_destroy_plan(fft);
  fftw_destroy_plan(ifft);
  fftw_free(aux);
  fftw_free(AUX);

  if (DBGFLG>1) {printf("main: proceed with hybrid extension\n"); fflush(NULL);}


  free(lsmask.lat);
  free(lsmask.lon);
  free2(lsmask.data, lsmask.nlon);
  free(freqs);

  return -1;
}
