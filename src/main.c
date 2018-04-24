#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
#include <netcdf.h>
#include <string.h>
#include "../lib/constants.h"
#include "../lib/dalloc.h"
#include "../lib/salloc.h"
#include "../lib/structs.h"
#include "save.h"
#include "damping.h"
#include "solveode.h"
#include "getdata.h"
#include "header.h"
#include "lsm.h"
#include "divergence.h"
#include "wavelength.h"
#include "hybrid.h"
#include "readtxt.h"
#include "../lib/macros.h"

int main(void) {
  if (DBGFLG > 2) {
    printf("main: get parameters from setup file\n");
    fflush(NULL);
  }

  // set expected amount of parameters AND paths (valLen) as well as amount of expected paths (pathLen)
  // Even if variables are set read-only (const) the compiler does not considers it a build time constant
  // and throws an error. Hence I defined them in a macro in header.h. They can be changed there.
  int valLen = VALLEN;
  int pathLen = PATHLEN;

  // iterator integers
  int nn, mm, ll;

  // allocate array of strings
  char** paths;
  paths = salloc2(paths, (size_t) pathLen, MAXCHARLEN);

  // initialize paths with null characters
  for (nn = 0; nn < pathLen; nn++) {
    for (mm = 0; mm < MAXCHARLEN; mm++) {
      paths[nn][mm] = '\0';
    }
  }

  // the length of the parameter array is exact. Increase if adding new (double) parameters
  double* params;
  params = dalloc(params, 17);
  readtxt(params, paths, valLen, pathLen);

  if (DBGFLG > 2) {
    printf("main: initialize variables\n");
    fflush(NULL);
  }

  // constants
  double r0, f0;

  // subsetting related variables
  int NLATMIN, NLATMAX, nlatmin, nlatmax, NLONMIN, NLONMAX, nlat5=NC_FILL_INT, slat5=NC_FILL_INT, hemflag;
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
  if ((int) params[10]%4==0){
    if ((int) params[10]%100!=0){
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

  // initialize variables with zeros
  for (nn=0; nn< (10248 + leap * 24); nn++){
    mld[nn] = taux[nn] = tauy[nn] = 0.0;
    time[nn] = 0;
  }

  // set basic time constructor
  tcon.tm_year = (int) params[10] - 1901;
  tcon.tm_mon = 11;
  tcon.tm_mday = 1;
  tcon.tm_hour = 1;
  tcon.tm_min = 0;
  tcon.tm_sec = 0;
  // tcon.tm_gmtoff = 0;
  tcon.tm_isdst = 0;

  // set reference time from time costructor
  time_t referenceTime;
  referenceTime = mktime(&tcon);

  // use reference time to set time axis
  for (nn=0; nn< ((365 + 62 + leap) * 24); nn++){
    time[nn] = (int) referenceTime + nn*3600;
  }

  if (DBGFLG > 1) {printf("main: check if domain is valid\n");fflush(NULL);}

  // if the domain is within the 5 degree band give error
  if (((int) params[11]>-5)&((int) params[12]<5)) {
    DOMAINERR
  }

  if (DBGFLG > 0) {printf("main: set land sea mask\n");fflush(NULL);}

  // allocate and inititate land sea mask
  dat2d lsmask = lsm(paths);

  if (DBGFLG>1) {printf("main: get subset\n"); fflush(NULL);}

  // find indices of subset in land sea mask
  NLATMIN = NLATMAX = NLONMIN = NLONMAX = 0;
  minimum = maximum = 180;

  // iterate over all latitudes
  for (nn=0; nn<lsmask.nlat; nn++){
    if (minimum>fabs(lsmask.lat[nn] - params[11])) {
      NLATMIN = nn;
      minimum = fabs(lsmask.lat[nn] - params[11]);
    }
    if (maximum>fabs(lsmask.lat[nn] - params[12])){
      NLATMAX = nn;
      maximum = fabs(lsmask.lat[nn] - params[12]);
    }
  }
  // correct if nearest neighbor is outside given interval
  if (NLATMAX<NLATMIN) {
    if (lsmask.lat[NLATMIN] < params[11]) NLATMIN -= 1;
    if (lsmask.lat[NLATMAX] > params[12]) NLATMAX += 1;
  } else {
    if (lsmask.lat[NLATMIN] < params[11]) NLATMIN += 1;
    if (lsmask.lat[NLATMAX] > params[12]) NLATMAX -= 1;
  }

  // the data requires positive longitudes within 0..360; correct if negative
  if (params[13] < 0) params[13] += 360;
  if (params[14] < 0) params[14] += 360;

  // iterate over all longitudes
  minimum = maximum = 360;
  for (nn=0; nn<lsmask.nlon; nn++){
    if (minimum>fabs(lsmask.lon[nn] - params[13])){
      NLONMIN = nn;
      minimum = fabs(lsmask.lon[nn] - params[13]);
    }
    if (maximum>fabs(lsmask.lon[nn] - params[14])){
      NLONMAX = nn;
      maximum = fabs(lsmask.lon[nn] - params[14]);
    }
  }
  // correct if nearest neighbor is outside given interval
  if (lsmask.lon[NLONMIN] < params[13]) NLONMIN +=1;
  if (lsmask.lon[NLONMAX] > params[14]) NLONMAX -=1;

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

  // set frequencies
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


  // be aware that the lat array may be sorted inversely
  if ((int) params[0] == 1) { // init nc files only if slab model is supposed to run.
    if (DBGFLG>2) {printf("main: init data files\n"); fflush(NULL);}
    initnc(params, paths, &lsmask, time, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, &nlat5, &slat5);
  } else if ((int) params[0] == 0) { // only load nlat5 and slat5; data should already be written.
    // set nlat5 and slat5
    for (nn=nlatmin; nn<nlatmax; nn++){
      if (dabs(lsmask.lat[nn]) < 5) break;
    }

    if (lsmask.lat[nn] > 0){
      memcpy(&nlat5, &nn, sizeof(int));
    } else {
      memcpy(&slat5, &nn, sizeof(int));
    }

    for (nn=nlatmax; nn>nlatmin; nn--){
      if (dabs(lsmask.lat[nn]) < 5) break;
    }

    if (lsmask.lat[nn] > 0){
      memcpy(&nlat5, &nn, sizeof(int));
    } else {
      memcpy(&slat5, &nn, sizeof(int));
    }

    // if nlat5 and slat5 did not receive a value throw an error
    // this is a heuristic check
    if ((nlat5==NC_FILL_INT)&(slat5==NC_FILL_INT)){
      GENERR
    }
  }

  // allocate and initialize the damping coefficient from Park et al 2009
  dat2d_2 damp = initdamping(&lsmask);

  // print settings of model before running
  if (DBGFLG>0) {
    printf("main: subset boundaries are set by:\n");
    printf("     (latmin, latmax, lonmin, lonmax):\n");
    printf("       %.2f; %.2f; %.2f; %.2f;\n",
           lsmask.lat[NLATMIN], lsmask.lat[NLATMAX], lsmask.lon[NLONMIN], lsmask.lon[NLONMAX]);
  }

  if ((int) params[0] == 1) {
    if (DBGFLG > 0) {
      printf("main: running slab model\n");
    }
    if (DBGFLG > 1) {
      printf("main: begin loop over points\n");
    }
    fflush(NULL);

    // calculate slab model on all points
    for (nn = nlatmin; nn <= nlatmax; nn++) {
      for (mm = NLONMIN; mm <= NLONMAX; mm++) {
        // check if point is on land; else leave fill value in file
        if (lsmask.data[nn][mm] == 0.0) {
          if (DBGFLG > 1) {
            printf("    (%d, %d)\n", nn, mm);
            fflush(NULL);
          }
          // load stress and mixed layer depth data
          getdata(params, paths, nn, mm, leap, time, mld, taux, tauy);
          // set damping parameter
          r0 = damping(&damp, nn, lsmask.lon[mm]);
          // set Coriolis frequency
          f0 = coriolis(lsmask.lat[nn]);
          // calculate u and v
          solve(fft, ifft, r0, f0, params[15], leap, taux, tauy, mld, freqs, aux, AUX);
          // set velocities
          for (ll = 0; ll < ((365 + 62 + leap) * 24); ll++) {
            uu[ll] = aux[ll][0];
            vv[ll] = aux[ll][1];
          }
          // write out mld, time, u and v
          savePoint(params, paths, &lsmask, uu, vv, mld, taux, tauy, time,
                    NLONMIN, mm, nlatmin, nlatmax, nn, leap, nlat5, slat5);
        }
      }
    }
  }

  if ((int) params[1] == 1) {
    if (DBGFLG > 0) {
      printf("main: run hybrid extension\n");
      fflush(NULL);
    }
    if (DBGFLG > 1) {
      printf("main: get vertical velocities and horizontal wavelengths\n");
      fflush(NULL);
    }

    // get mid-point divergences and horizontal wavelengths on both hemispheres
    // parameters are set according to domain setup
    // possible setups:
    //  - minimum and maximum latitude on northern hemisphere
    //  - minimum and maximum latitude on southern hemisphere
    //  - the minimum latitude is in the forbidden 5 degree band
    //  - the maximum latitude is in the forbidden 5 degree band
    //  - full case with two valid hemispheres
    // this setups are possible with descending or ascending latitude array

    if (NLATMIN == nlatmin) { // ascending latitude array
      if (((int) params[11] > 0) & ((int) params[12] > 0)) { // both on northern hemisphere

        if ((int) params[2] == 1) divergence(&lsmask, params, paths, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, 0);
        if ((int) params[3] == 1) wavelength(&lsmask, params, paths, time, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, 0);

      } else if (((int) params[11] < 0) & ((int) params[12] < 0)) { // both on southern hemisphere

        if ((int) params[2] == 1) divergence(&lsmask, params, paths, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, 1);
        if ((int) params[3] == 1) wavelength(&lsmask, params, paths, time, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, 1);

      } else if ((int) params[11] > -5) { // latitude minimum in forbidden band

        if ((int) params[2] == 1) divergence(&lsmask, params, paths, NLONMIN, NLONMAX + 1, nlat5, nlatmax + 1, leap, 0);
        if ((int) params[3] == 1) wavelength(&lsmask, params, paths, time, NLONMIN, NLONMAX + 1, nlat5, nlatmax + 1, leap, 0);

      } else if ((int) params[12] < 5) { // latitude maximum in forbidden band

        if ((int) params[2] == 1) divergence(&lsmask, params, paths, NLONMIN, NLONMAX + 1, nlatmin, slat5 + 1, leap, 1);
        if ((int) params[3] == 1) wavelength(&lsmask, params, paths, time, NLONMIN, NLONMAX + 1, nlatmin, slat5 + 1, leap, 1);

      } else { // general case with two valid hemispheres

        if ((int) params[2] == 1) divergence(&lsmask, params, paths, NLONMIN, NLONMAX + 1, nlat5, nlatmax + 1, leap, 0);
        if ((int) params[3] == 1) wavelength(&lsmask, params, paths, time, NLONMIN, NLONMAX + 1, nlat5, nlatmax + 1, leap, 0);

        if ((int) params[2] == 1) divergence(&lsmask, params, paths, NLONMIN, NLONMAX + 1, nlatmin, slat5 + 1, leap, 1);
        if ((int) params[3] == 1) wavelength(&lsmask, params, paths, time, NLONMIN, NLONMAX + 1, nlatmin, slat5 + 1, leap, 1);

      }
    } else if (NLATMAX == nlatmin) { // descending latitude array
      if (((int) params[11] > 0) & ((int) params[12] > 0)) { // both on northern hemisphere

        /* for (nn=0; nn < lsmask.nlat; nn++) {
          printf("nn, lat[nn]: %d\t%f\n", nn, lsmask.lat[nn]);fflush(NULL);
        } */

        if ((int) params[2] == 1) divergence(&lsmask, params, paths, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, 0);
        if ((int) params[3] == 1) wavelength(&lsmask, params, paths, time, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, 0);

      } else if (((int) params[11] < 0) & ((int) params[12] < 0)) { // both on southern hemisphere

        if ((int) params[2] == 1) divergence(&lsmask, params, paths, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, 1);
        if ((int) params[3] == 1) wavelength(&lsmask, params, paths, time, NLONMIN, NLONMAX + 1, nlatmin, nlatmax + 1, leap, 1);

      } else if ((int) params[11] > -5) { // latitude minimum in forbidden band

        if ((int) params[2] == 1) divergence(&lsmask, params, paths, NLONMIN, NLONMAX + 1, nlatmin, nlat5 + 1, leap, 0);
        if ((int) params[3] == 1) wavelength(&lsmask, params, paths, time, NLONMIN, NLONMAX + 1, nlatmin, nlat5 + 1, leap, 0);

      } else if ((int) params[12] < 5) { // latitude maximum in forbidden band

        if ((int) params[2] == 1) divergence(&lsmask, params, paths, NLONMIN, NLONMAX + 1, slat5, nlatmax + 1, leap, 1);
        if ((int) params[3] == 1) wavelength(&lsmask, params, paths, time, NLONMIN, NLONMAX + 1, slat5, nlatmax + 1, leap, 1);

      } else { // general case with two valid hemispheres

        if ((int) params[2] == 1) divergence(&lsmask, params, paths, NLONMIN, NLONMAX + 1, nlatmin, nlat5 + 1, leap, 0);
        if ((int) params[3] == 1) wavelength(&lsmask, params, paths, time, NLONMIN, NLONMAX + 1, nlatmin, nlat5 + 1, leap, 0);

        if ((int) params[2] == 1) divergence(&lsmask, params, paths, NLONMIN, NLONMAX + 1, slat5, nlatmax + 1, leap, 1);
        if ((int) params[3] == 1) wavelength(&lsmask, params, paths, time, NLONMIN, NLONMAX + 1, slat5, nlatmax + 1, leap, 1);

      }
    }

    // set structs for time series of N, w, horiz. wavelength lh and radiated energy Eout
    dat1d lh;
    dat1d ww;
    dat1d NN;
    dat1d Eout;
    lh.ntime = (365 + leap) * 24;
    NN.ntime = ww.ntime = Eout.ntime = lh.ntime;
    lh.time = dalloc(lh.time, (size_t) lh.ntime);
    for (nn = 0; nn < lh.ntime; nn++) lh.time[nn] = (double) time[nn + 31 * 24];
    NN.time = ww.time = Eout.time = lh.time;

    // allocate data
    lh.data = dalloc(lh.data, (size_t) lh.ntime);
    ww.data = dalloc(lh.data, (size_t) lh.ntime);
    NN.data = dalloc(lh.data, (size_t) lh.ntime);
    Eout.data = dalloc(lh.data, (size_t) lh.ntime);

    if (DBGFLG > 1) {
      printf("main: begin loop over points\n");
      fflush(NULL);
    }

    // calculate hybrid model on all points
    for (nn = nlatmin + 1; nn < nlatmax; nn++) {
      for (mm = NLONMIN + 1; mm < NLONMAX; mm++) {
        // check if point has neighbors, if not skip
        if (lsmask.data[nn - 1][mm] +
            lsmask.data[nn + 1][mm] +
            lsmask.data[nn][mm - 1] +
            lsmask.data[nn][mm + 1] +
            lsmask.data[nn][mm] == 0.0) {
          if (DBGFLG > 1) {
            printf("    (%d, %d)\n", nn, mm);
            fflush(NULL);
          }

          // get data
          getdataHybrid(params, paths, &lsmask, &lh, &ww, &NN, nn, nlatmin, mm, NLONMIN, leap, nlat5, slat5);

          if (lsmask.lat[nn]<-5){
            // this is a debugging test
            printf("southern hemisphere, data test");fflush(NULL);
          }

          // set Coriolis frequency
          f0 = coriolis(lsmask.lat[nn]);

          // get energy flux for point
          hybrid(params, &lh, &ww, &NN, &Eout, hfreqs, f0, HAUX, haux, hfft, hifft, leap);

          // save data to nc file
          savePointHybrid(params, paths, &lsmask, &Eout, nn, nlatmin, nlatmax, mm, NLONMIN, leap, nlat5, slat5);
        }
      }
    }

    // free up memory from data used in hybrid extension
    free(lh.data);
    free(lh.time); // frees time axes of all data structs
    free(NN.data);
    free(ww.data);
    free(Eout.data);
  }

  if (DBGFLG>2) {printf("main: clean up memory\n");fflush(NULL);}

  // free up memory for clean exit
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
  free(params);
  sfree2(paths, (size_t) pathLen);

  // if you made it here exit with success
  return(EXIT_SUCCESS);
}