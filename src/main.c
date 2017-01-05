#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
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
#include "wavelength.h"
#include "hybrid.h"

int main(void) {
  if (DBGFLG > 2) {
    printf("main: initialize variables\n");
    fflush(NULL);
  }
  // iterator integers
  int nn, mm, ll;

  // constants
  double rho0 = RHO, r0, f0;

  // subsetting related variables
  int NLATMIN, NLATMAX, nlatmin, nlatmax, NLONMIN, NLONMAX, nlat5=0, slat5=0, hemflag;
  double minimum, maximum;

  // time loop related variables
  int leap;
  struct tm tcon;

  // FFT related variables
  double freqs[DFFT_LEN], *hfreqs;
  fftw_complex *aux, *AUX, *haux, *HAUX;
  fftw_plan fft, ifft, hfft, hifft;

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
  double taux[10248 + leap * 24], tauy[10248 + leap * 24],
          mld[10248 + leap * 24],
          uu[10248 + leap * 24], vv[10248 + leap * 24];
  int time[10248 + leap * 24];

  for (nn=0; nn< (10248 + leap * 24); nn++){
    mld[nn] = taux[nn] = tauy[nn] = 0.0;
    time[nn] = 0;
  }

  // set basic time constructor
  tcon.tm_year = YEAR - 1901;
  tcon.tm_mon = 11;
  tcon.tm_mday = 1;
  tcon.tm_hour = 1;
  tcon.tm_min = 0;
  tcon.tm_sec = 0;
  // tcon.tm_gmtoff = 0;
  tcon.tm_isdst = 0;

  time_t referenceTime;
  referenceTime = mktime(&tcon);
  for (nn=0; nn< ((365 + 62 + leap) * 24); nn++){
    time[nn] = (int) referenceTime + nn*3600;

  }

  if (DBGFLG > 1) {printf("main: check if domain is valid\n");fflush(NULL);}

  if ((LATMIN>-5)&(LATMAX<5)) {
    DOMAINERR
  }

  if (DBGFLG > 0) {printf("main: set land sea mask\n");fflush(NULL);}
  dat2d lsmask = lsm();

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
  haux = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (365+leap)*24);
  HAUX = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (365+leap)*24);
  hfreqs = dalloc(hfreqs, (size_t) ((365 + leap) * 24));

  for (nn=0; nn<DFFT_LEN/2; nn++) {
    freqs[nn] = 2 * PI / 3600.0 * nn / DFFT_LEN;
  }
  for (nn=DFFT_LEN/2; nn<DFFT_LEN; nn++) {
    freqs[nn] = 2 * PI / 3600.0 * (nn - DFFT_LEN) / DFFT_LEN;
  }

  for (nn=0; nn<((365 + leap) * 24) / 2; nn++) {
    hfreqs[nn] = 2 * PI / 3600.0 * nn / ((365 + leap) * 24);
  }
  for (nn=((365 + leap) * 24) / 2; nn<((365 + leap) * 24); nn++) {
    hfreqs[nn] = 2 * PI / 3600.0 * (nn - ((365 + leap) * 24)) / ((365 + leap) * 24);
  }

  // define transforms
  fft = fftw_plan_dft_1d(DFFT_LEN, aux, AUX, FFTW_FORWARD, FFTW_ESTIMATE);
  ifft = fftw_plan_dft_1d(DFFT_LEN, AUX, aux, FFTW_BACKWARD, FFTW_ESTIMATE);

  hfft = fftw_plan_dft_1d(((365 + leap) * 24), haux, HAUX, FFTW_FORWARD, FFTW_ESTIMATE);
  hifft = fftw_plan_dft_1d(((365 + leap) * 24), HAUX, haux, FFTW_BACKWARD, FFTW_ESTIMATE);

  if (DBGFLG>2) {printf("main: init data files\n"); fflush(NULL);}

  // be aware that the lat array may be sorted inversely
  initnc(&lsmask, time, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, nlat5, slat5);
  dat2d_2 damp = initdamping(&lsmask);

  if (DBGFLG>0) {
    printf("main: subset boundaries are set by:\n");
    printf("     (latmin, latmax, lonmin, lonmax):\n");
    printf("       %.2f; %.2f; %.2f; %.2f;\n",
           lsmask.lat[NLATMIN], lsmask.lat[NLATMAX], lsmask.lon[NLONMIN], lsmask.lon[NLONMAX]);
    printf("main: running slab model\n");
  }
  if (DBGFLG>1) {
    printf("main: begin loop over points\n");
  }
  fflush(NULL);

  // calculate slab model on all points
  for (nn=nlatmin; nn<=nlatmax; nn++){
    for (mm=NLONMIN; mm<=NLONMAX; mm++){
      // check if point is on land; else leave fill value in file
      if (lsmask.data[nn][mm]==0.0) {
        if (DBGFLG > 1) {
          printf("    (%d, %d)\n", nn, mm);
          fflush(NULL);
        }
        // load stress and mixed layer depth data
        getdata(nn, mm, leap, time, mld, taux, tauy);
        // set damping parameter
        r0 = damping(&damp, nn, lsmask.lon[mm]);
        // set Coriolis frequency
        f0 = coriolis(lsmask.lat[nn]);
        // calculate u and v
        solve(fft, ifft, r0, f0, rho0, leap, taux, tauy, mld, freqs, aux, AUX);
        // set velocities
        for (ll = 0; ll < ((365 + 62 + leap) * 24); ll++) {
          uu[ll] = aux[ll][0];
          vv[ll] = aux[ll][1];
        }
        // write out mld, time, u and v
        savePoint(&lsmask, uu, vv, mld, taux, tauy, time, NLONMIN, mm, nlatmin, nn, leap, nlat5, slat5);
      }
    }
  }

  if (DBGFLG>0) {printf("main: proceed with hybrid extension\n"); fflush(NULL);}

  if (DBGFLG>1) {printf("main: get vertical velocities and horizotal wavelengths\n");fflush(NULL);}

  // get mid-point divergences on both hemispheres

  if ((LATMIN > 0) & (LATMAX > 0)) { // both on northern hemisphere

    divergence(&lsmask, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, 0);
    wavelength(&lsmask, time, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, 0);

  } else if ((LATMIN < 0) & (LATMAX < 0)) { // both on southern hemisphere

    divergence(&lsmask, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, 1);
    wavelength(&lsmask, time, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, 1);

  }

  if (NLATMIN == nlatmin) {
    if (LATMIN > -5) { // latitude minimum in forbidden band

      divergence(&lsmask, NLONMIN, NLONMAX + 1, nlat5, nlatmax + 1, leap, 0);
      wavelength(&lsmask, time, NLONMIN, NLONMAX + 1, nlat5, nlatmax + 1, leap, 0);

    } else if (LATMAX < 5) { // latitude maximum in forbidden band

      divergence(&lsmask, NLONMIN, NLONMAX + 1, nlatmin, slat5 + 1, leap, 1);
      wavelength(&lsmask, time, NLONMIN, NLONMAX + 1, nlatmin, slat5 + 1, leap, 1);

    } else { // general case with two valid hemispheres

      divergence(&lsmask, NLONMIN, NLONMAX + 1, nlat5, nlatmax + 1, leap, 0);
      wavelength(&lsmask, time, NLONMIN, NLONMAX + 1, nlat5, nlatmax + 1, leap, 0);

      divergence(&lsmask, NLONMIN, NLONMAX + 1, nlatmin, slat5 + 1, leap, 1);
      wavelength(&lsmask, time, NLONMIN, NLONMAX + 1, nlatmin, slat5 + 1, leap, 1);

    }
  } else if (NLATMAX == nlatmin) {
    if (LATMIN > -5) { // latitude minimum in forbidden band

      divergence(&lsmask, NLONMIN, NLONMAX + 1, nlatmin, nlat5 + 1, leap, 0);
      wavelength(&lsmask, time, NLONMIN, NLONMAX + 1, nlatmin, nlat5 + 1, leap, 0);

    } else if (LATMAX < 5) { // latitude maximum in forbidden band

      divergence(&lsmask, NLONMIN, NLONMAX + 1, slat5, nlatmax + 1, leap, 1);
      wavelength(&lsmask, time, NLONMIN, NLONMAX + 1, slat5, nlatmax + 1, leap, 1);

    } else { // general case with two valid hemispheres

      divergence(&lsmask, NLONMIN, NLONMAX + 1, nlatmin, nlat5 + 1, leap, 0);
      wavelength(&lsmask, time, NLONMIN, NLONMAX + 1, nlatmin, nlat5 + 1, leap, 0);

      divergence(&lsmask, NLONMIN, NLONMAX + 1, slat5, nlatmax + 1, leap, 1);
      wavelength(&lsmask, time, NLONMIN, NLONMAX + 1, slat5, nlatmax + 1, leap, 1);

    }
  }



  // set struct for data time series
  dat1d lh;
  dat1d ww;
  dat1d NN;
  dat1d Eout;
  lh.ntime = (365+leap)*24;
  NN.ntime = ww.ntime = Eout.ntime = lh.ntime;
  lh.time = dalloc(lh.time, (size_t) lh.ntime);
  for (nn=0; nn < lh.ntime; nn++) lh.time[nn] = (double) time[nn + 31*24];
  NN.time = ww.time = Eout.time = lh.time;

  lh.data = dalloc(lh.data, (size_t) lh.ntime);
  ww.data = dalloc(lh.data, (size_t) lh.ntime);
  NN.data = dalloc(lh.data, (size_t) lh.ntime);
  Eout.data = dalloc(lh.data, (size_t) lh.ntime);

  if (DBGFLG>1) {printf("main: begin loop over points\n");fflush(NULL);}
  int nt;
  // calculate hybrid model on all points
  for (nn=nlatmin +1 ; nn<nlatmax; nn++){
    for (mm=NLONMIN + 1; mm<NLONMAX; mm++){

      if (lsmask.data[nn-1][mm] +
          lsmask.data[nn+1][mm] +
          lsmask.data[nn][mm-1] +
          lsmask.data[nn][mm+1] +
          lsmask.data[nn][mm] == 0.0) {
        if (DBGFLG>1) { printf("    (%d, %d)\n", nn, mm); fflush(NULL);}


        // get data
        getdataHybrid(&lsmask, &lh, &ww, &NN, nn, nlatmin, mm, NLONMIN, leap, nlat5, slat5);

        // set Coriolis frequency
        f0 = coriolis(lsmask.lat[nn]);

        // get energy flux for point
        hybrid(&lh, &ww, &NN, &Eout, hfreqs, f0, HAUX, haux, hfft, hifft, leap);

        // save data to nc file
        savePointHybrid(&lsmask, &Eout, nn, nlatmin, mm, NLONMIN, leap, nlat5, slat5);
      }
    }
  }

  if (DBGFLG>2) {printf("main: clean up memory\n");fflush(NULL);}

  fftw_destroy_plan(fft);
  fftw_destroy_plan(ifft);
  fftw_destroy_plan(hfft);
  fftw_destroy_plan(hifft);
  fftw_free(aux);
  fftw_free(AUX);
  free(lsmask.lat);
  free(lsmask.lon);
  dfree2(lsmask.data, (size_t) lsmask.nlat);
  free(damp.y1); // note tha damp.xx is free'd at lsmask.lat
  free(damp.y2);
  free(lh.data);
  free(lh.time); // frees time axes of all data structs
  free(NN.data);
  free(ww.data);
  free(Eout.data);

  return(EXIT_SUCCESS);
}