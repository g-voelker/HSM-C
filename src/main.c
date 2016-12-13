#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
//#include "constants.h"
//#include "macros.h"
#include "../lib/constants.h"
#include "../lib/dalloc.h"
#include "../lib/structs.h"
#include "save.h"
#include "damping.h"
#include "solveode.h"
#include "getdata.h"
#include "header.h"
#include "lsm.h"
#include "input.h"
#include "divergence.h"

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
  if (DBGFLG > 2) {
    printf("main: initialize variables\n");
    fflush(NULL);
  }
  // iterator integers
  int nn, mm, ll;

  // constants
  double rho0 = RHO, r0, f0;

  // land sea mask related variables
  char lsmfile[] = LSMPATH;

  // subsetting related variables
  int NLATMIN, NLATMAX, nlatmin, nlatmax, NLONMIN, NLONMAX;
  double minimum, maximum;

  // time loop related variables
  int leap;
  struct tm tcon;

  // FFT related variables
  double freqs[DFFT_LEN];
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
          mld[8760 + leap * 24],
          uu[8760 + leap * 24], vv[8760 + leap * 24];
  int time[8760 + leap * 24];

  for (nn=0; nn< (8760 + leap * 24); nn++){
    mld[nn] = taux[nn] = tauy[nn] = 0.0;
    time[nn] = 0;
  }

  // set basic time constructor
  tcon.tm_year = YEAR - 1900;
  tcon.tm_mon = 0;
  tcon.tm_mday = 1;
  tcon.tm_hour = 1;
  tcon.tm_min = 0;
  tcon.tm_sec = 0;
  tcon.tm_gmtoff = 0;
  tcon.tm_isdst = 0;

  time_t referenceTime;
  referenceTime = mktime(&tcon);
  for (nn=0; nn< ((365 + leap) * 24); nn++){
    time[nn] = (int) referenceTime + nn*3600;

  }

  if (DBGFLG > 1) {printf("main: set land sea mask\n");fflush(NULL);}
  dat2d lsmask = lsm(lsmfile);

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
  // correct if nearest neighbor is outside given interval
  if (NLATMAX<NLATMIN) {
    if (lsmask.lat[NLATMIN] < LATMIN) NLATMIN -= 1;
    if (lsmask.lat[NLATMAX] > LATMAX) NLATMAX += 1;
  } else {
    if (lsmask.lat[NLATMIN] < LATMIN) NLATMIN += 1;
    if (lsmask.lat[NLATMAX] > LATMAX) NLATMAX -= 1;
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
  // correct if nearest neighbor is outside given interval
  if (lsmask.lon[NLONMIN] < LONMIN) NLONMIN +=1;
  if (lsmask.lon[NLONMAX] > LONMAX) NLONMAX -=1;

  // lat array may be sorted reversely
  nlatmin = (((NLATMIN)<(NLATMAX))?(NLATMIN):(NLATMAX));
  nlatmax = (((NLATMIN)>(NLATMAX))?(NLATMIN):(NLATMAX));

  if (DBGFLG>2) {printf("main: set fftw plans\n"); fflush(NULL);}
  // prepare FFTs
  // set frequencies
  aux = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * DFFT_LEN);
  AUX = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * DFFT_LEN);

  for (nn=0; nn<DFFT_LEN/2; nn++) {
    freqs[nn] = 2 * PI / 3600.0 * nn / DFFT_LEN;
  }
  for (nn=DFFT_LEN/2; nn<DFFT_LEN; nn++) {
    freqs[nn] = 2 * PI / 3600.0 * (nn - DFFT_LEN) / DFFT_LEN;
  }

  // define transforms
  fft = fftw_plan_dft_1d(DFFT_LEN, aux, AUX, FFTW_FORWARD, FFTW_ESTIMATE);
  ifft = fftw_plan_dft_1d(DFFT_LEN, AUX, aux, FFTW_BACKWARD, FFTW_ESTIMATE);

  if (DBGFLG>2) {printf("main: init data files\n"); fflush(NULL);}

  // be aware that the lat array may be sorted inversely
  initnc(&lsmask, time, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap);
  dat2d_2 damp = initdamping(&lsmask);

  if (DBGFLG>1) {
    printf("main: subset boundaries are set by:\n");
    printf("     (latmin, latmax, lonmin, lonmax):\n");
    printf("       %.2f; %.2f; %.2f; %.2f;\n",
           lsmask.lat[NLATMIN], lsmask.lat[NLATMAX], lsmask.lon[NLONMIN], lsmask.lon[NLONMAX]);
    printf("main: begin loop over points\n");
    fflush(NULL);
  }

  fflush(NULL);

  // calculate slab model on all points
  for (nn=nlatmin; nn<=nlatmax; nn++){
    for (mm=NLONMIN; mm<=NLONMAX; mm++){
      if (DBGFLG>1) { printf("    (%d, %d)\n", nn, mm); fflush(NULL);}
      // load stress and mixed layer depth data
      getdata(nn, mm, leap, time, mld, taux, tauy);
      // set damping parameter
      r0 = damping(&damp, nn, lsmask.lon[mm]);
      // set Coriolis frequency
      f0 = coriolis(lsmask.lat[nn]);
      // calculate u and v
      solve(fft, ifft, r0, f0, rho0, leap, taux, tauy, mld, freqs, aux, AUX);
      // set velocities
      for (ll=0; ll < ((365 + leap) * 24); ll++){
        uu[ll] = aux[ll][0];
        vv[ll] = aux[ll][1];
      }
      // write out mld, time, u and v
      savePoint(uu, vv, mld, taux, tauy, time, NLONMIN, mm, nlatmin, nn, leap);
    }
  }

  if (DBGFLG>1) {printf("main: proceed with hybrid extension\n"); fflush(NULL);}

  // get mid-point divergences

  divergence(&lsmask, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap);

  // get energy fluxes


  if (DBGFLG>2) {printf("main: clean up memory\n");fflush(NULL);}

  fftw_destroy_plan(fft);
  fftw_destroy_plan(ifft);
  fftw_free(aux);
  fftw_free(AUX);
  free(lsmask.lat);
  free(lsmask.lon);
  dfree2(lsmask.data, (size_t) lsmask.nlon);

  return -1;
}
