//
// Created by georg on 14/12/16.
//

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <math.h>
#include <string.h>
#include "../lib/structs.h"
#include "../lib/dalloc.h"
#include "../lib/stats.h"
#include "input.h"
#include "header.h"
#include "save.h"
#include "divergence.h"

void loadw(dat3d *ww, dat2d *lsmask, int nlatmin, int nlon, int ndlon, int nmonth, int leap){
  // load slice of vertical velocity from data files
  int retval, ncID, varID, nn, mm;
  size_t days[12] = {31, 28 + (size_t) leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  size_t start[3] = {0, 0, 0}, count[3]={0, 1, 1}; // time, lats, lons
  char filepath[MAXCHARLEN];

  // set file to load from
  if (nmonth==0) {
    sprintf(filepath, AUXPATH, 1);
    count[0] = days[0] * 24;
  } else if (nmonth==13) {
    sprintf(filepath, AUXPATH, 12);
    count[0] = days[11] * 24;
  } else {
    sprintf(filepath, OUTPATH, nmonth);
    count[0] = days[nmonth-1] * 24;
  }

  // count indicees
  start[0] = 0;

  if (DBGFLG>2) {printf("  loadw: get first slice of vertical velocity\n");fflush(NULL);}

  // open file
  if ((retval = nc_open(filepath, NC_NOWRITE, &ncID))) ERR(retval);

//  int dimid[3] = {0,0,0};
//  size_t dimsize[3] = {0,0,0};
//
//  if ((retval = nc_inq_dimid(ncID, LATS, &dimid[0]))) ERR(retval);
//  if ((retval = nc_inq_dimid(ncID, LONS, &dimid[1]))) ERR(retval);
//  if ((retval = nc_inq_dimid(ncID, TIME, &dimid[2]))) ERR(retval);
//
//  if ((retval = nc_inq_dimlen(ncID, dimid[0], &dimsize[0]))) ERR(retval);
//  if ((retval = nc_inq_dimlen(ncID, dimid[1], &dimsize[1]))) ERR(retval);
//  if ((retval = nc_inq_dimlen(ncID, dimid[2], &dimsize[2]))) ERR(retval);
//  printf("dimension sizes: (%d, %d, %d)\n", (int) dimsize[0], (int) dimsize[1], (int) dimsize[2]);fflush(NULL);

  // get varID of vertical velocity
  if ((retval = nc_inq_varid(ncID, ZVEL, &varID))) ERR(retval);

  // read slices
  for (nn=nlatmin; nn<nlatmin + ww->nlat; nn++){
    start[1] = (size_t) (nn - nlatmin);
    for (mm=nlon; mm < nlon + 2 * ndlon; mm++) {
      // read slice
      if (mm < 0) {
        start[2] = (size_t) (lsmask->nlon + mm);
      } else if (mm > lsmask->nlon) {
        start[2] = (size_t) (mm - lsmask->nlon);
      }

      if ((retval = nc_get_vara_double(ncID, varID, start, count, &(ww->data)[nn-nlatmin][mm - nlon][0])))

      ERR(retval);
    }
  }

  // close nc file
  if ((retval = nc_close(ncID))) ERR(retval);
}

void advancew(dat3d *ww, dat2d *lsmask, int nlatmin, int nlon, int ndlon, int nmonth, int leap){
  // advance by one index
  // load slice of vertical velocity from data files
  int retval, ncID, varID, nn, mm, nt;
  size_t days[12] = {31, 28 + (size_t) leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  size_t start[3] = {0, 0, 0}, count[3]={0, 1, 1}; // time, lats, lons
  char filepath[MAXCHARLEN];

  if (DBGFLG>2) {printf("  advancew: shifting vertical velocity\n");fflush(NULL);}
  // shift data by one
  for (nn=0; nn < ww->nlat; nn++){
    for (mm=0; mm < 2 * ndlon - 1; mm++){
      for (nt=0; nt < ww->nt; nt++){
        memcpy(&ww->data[nn][mm][nt], &ww->data[nn][mm+1][nt] , sizeof(double));
      }
    }
  }

  if (DBGFLG>2) {printf("  advancew: advance vertical velocity field by one longitude\n");fflush(NULL);}
  // load one slice of data
  // set file to load from
  if (nmonth==0) {
    sprintf(filepath, AUXPATH, 1);
    count[0] = days[0] * 24;
  } else if (nmonth==13) {
    sprintf(filepath, AUXPATH, 12);
    count[0] = days[11] * 24;
  } else {
    sprintf(filepath, OUTPATH, nmonth);
    count[0] = days[nmonth-1] * 24;
  }

  // open file
  if ((retval = nc_open(filepath, NC_NOWRITE, &ncID))) ERR(retval);

  // get varID of vertical velocity
  if ((retval = nc_inq_varid(ncID, ZVEL, &varID))) ERR(retval);

  // read arrays
  for (nn=nlatmin; nn<nlatmin + ww->nlat; nn++){
    mm = nlon + 2 * ndlon - 1;
    // read slice
    if (mm < 0) {
      start[2] = (size_t) (lsmask->nlon + mm);
    } else if (mm > lsmask -> nlon) {
      start[2] = (size_t) (mm - lsmask->nlon);
    }
    if ((retval = nc_get_vara_double(ncID, varID, start, count, &(ww->data)[nn-nlatmin][mm - nlon][0])))
    ERR(retval);
  }

  // close nc file
  if ((retval = nc_close(ncID))) ERR(retval);
}

void autocorr(dat2d *lsmask, dat3d *ww, double ***distances, double **lh,
              int nlon, int nlonmin, int nlonmax, int nlatmin,
              int leap, int ndlat, int ndlon, int nmonth, int edgeflag) {
  // declare used variables
  size_t days[14] = {31, 31, 28 + (size_t) leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31};
  int nn, mm, nx, ny, nt, index, iterbounds[4], iternlon[2], itershift;
  double *distance, *corr, *norm, avg;

  distance = dalloc(distance, CORRLEN);
  corr = dalloc(corr, CORRLEN);
  norm = dalloc(norm, CORRLEN);

  for (index=0; index<CORRLEN; index++) distance[index] = CORRMIN + (CORRMAX - CORRMIN)/(CORRLEN - 1) * index;

  if ((edgeflag==0)||((nlon >= nlonmin)&&(nlon <= nlonmax - 2 * ndlon))) {
    iternlon[0] = nlon + ndlon;
    iternlon[1] = nlon + ndlon + 1;
  } else if ((edgeflag==1)&&(nlon == nlonmin)) {
    iternlon[0] = nlonmin;
    iternlon[1] = nlonmin + ndlon + 1;
  } else if ((edgeflag==1)&&(nlon == nlonmax - 2*ndlon)) {
    iternlon[0] = nlon + ndlon;
    iternlon[1] = nlonmax;
  } else {
    GENERR; fflush(NULL);
  }

  for (mm=iternlon[0]; mm<iternlon[1]; mm++) {

    for (nn = nlatmin; nn < nlatmin + ww->nlat; nn++) {
      // initialize wavelength with zero
      lh[nn - nlatmin][mm - nlonmin] = 0.0;
      // boundaries for iteration over latitudes
      iterbounds[0] = (int) fmax(nn - ndlat, nlatmin);
      iterbounds[1] = (int) fmin(nn + ndlat, nlatmin + ww->nlat);
      iterbounds[2] = (int) fmax(mm - 2*ndlon, nlon);
      iterbounds[3] = (int) fmin(mm + ndlon, nlon + ww->nlon);

      // iterate over all time steps
      for (nt = 0; nt < days[nmonth]; nt++) {
        // get average vertical velocity
        avg = wavg2(lsmask, ww, nt, iterbounds[2], iterbounds[0]);
        // set correlation array and norms to zero
        for (index = 0; index < CORRLEN; index++) corr[index] = norm[index] = 0.0;
        // iterate over lats and lons and add element to corresponding element in correlation array
        // note that we neglect any factor constant at the point since we evaluate maxima only
        for (ny = iterbounds[0]; ny < iterbounds[1]; ny++) {
          for (nx = iterbounds[2]; nx < iterbounds[3]; nx++) {
            if (ww->data[nn - nlatmin][mm - nlon][nt]!=NC_FILL_DOUBLE) {
              // check distance and get index
              itershift = abs(mm - nlon - ndlon);
              index = (int) (
                      (distances[nn - nlatmin][ny - iterbounds[0]][nx - nlon + itershift] - CORRMIN) / (CORRMAX - CORRMIN) *
                      CORRLEN + 0.5);
              // add deviation to corr array and increase norm by 1
              if ((index < CORRLEN) & (index >= 0) & (ww->data[ny - nlatmin][nx - nlon][nt]!=NC_FILL_DOUBLE)) {
                corr[index] +=
                        (ww->data[ny - nlatmin][nx - nlon][nt] - avg) * (ww->data[nn - nlatmin][mm - nlon][nt] - avg);
                norm[index]++;
              }
            }
          }
        }
        // normalize correlation with number of points in bin
        for (index = 0; index < CORRLEN; index++) corr[index] /= norm[index];
        // get wavelegth from correlation array of point and add result to lh[nn-nlatmin][nlon-nlonmin]
        lh[nn - nlatmin][nlon - nlonmin] += dxmax(distance, corr, CORRLEN);
      }
      // normalize actual point in lh
      lh[nn - nlatmin][nlon - nlonmin] /= days[nmonth] * 24;
    }
  }

  free(distance);
  free(corr);
  free(norm);
}

void wavelength(dat2d *lsmask, int *time, int nxmin, int nxmax, int nymin, int nymax, int leap){
  size_t days[14] = {31, 31, 28 + (size_t) leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31};
  dat3d ww;
  size_t ntime;
  int nmonth, nlat, nlon, ndlon, nn, nx, ny, mm, itermax, edgeflag, ndlat, iterbounds[2];
  double latmax, dx;
  double ***distances, ***lh;

  // prepare constants
  nlat = nymax - nymin; // number of points in latuitude of w-slice
  latmax = fmax(fabs(lsmask->lat[nymin]), fabs(lsmask->lat[nymax])); // maximum latitude
  dx = DEG2RAD(lsmask->lon[1] - lsmask->lon[0])*EARTHRADIUS*cos(DEG2RAD(latmax)); // highest x-resolution in domain
  ndlon = (int) (CORRMAX/dx) + 2; // minimum number of points needed to include length from CORRMAX
  ndlat = (int) (CORRMAX / DEG2RAD(fabs(lsmask->lat[1] - lsmask->lat[0])) / EARTHRADIUS )  + 2;

  // if domain is smaller than needed give error and abort calculation
  // if ( ndlon > nxmax - nxmin) VALERR(2);
  nlon = 2 * ndlon; // number of points in lon of w-slice

  // check for repeating domain ('around the globe')
  edgeflag = 1;
  itermax = nxmax - 2 * ndlon;
  if ((nxmax - nxmin) == lsmask->nlon) {
    itermax = lsmask->nlon;
    edgeflag = 0;
  };

  if (DBGFLG>2) {printf("  wavelength: calculate distances\n");fflush(NULL);}
  // set up distances arrays
  iterbounds[0] = (int) fmin(2*ndlat, nymax - nymin);
  iterbounds[1] = (int) fmin(2*ndlon, nxmax - nxmin);
  distances = dalloc3(distances, (size_t) (nymax-nymin), (size_t) iterbounds[0], (size_t) iterbounds[1]);
  for (nn=nymin; nn<nymax; nn++){
    for (ny=0; ny<iterbounds[0]; ny++){
      for (nx=0; nx<iterbounds[1]; nx++){
        distances[nn - nymin][ny][nx] = dist(lsmask->lon[nx + nxmin], lsmask->lon[ndlon + nxmin],
                                             lsmask->lat[ny + nymin], lsmask->lat[nn]);
      }
    }
  }

  if (DBGFLG>2) {printf("  wavelength: get auto-correlation for all months\n");fflush(NULL);}
  // iterate over months
  lh = dalloc3(lh, 14, (size_t) nlat, (size_t) (nxmax-nxmin));
  for (nmonth=0; nmonth<14; nmonth++) {
    ntime = days[nmonth] * 24;
    ww.data = dalloc3(ww.data, (size_t) nlat, (size_t) nlon, ntime);
    ww.nlat = nlat;
    ww.nlon = nlon;
    ww.nt = (int) ntime;
    // iterate over longitudes (+ handle edges of domains)
    for (mm=nxmin; mm < itermax; mm++){
      // load vertical velocity band
      if (mm==nxmin){
        loadw(&ww, lsmask, nymin, mm, ndlon, nmonth, leap);
      } else {
        advancew(&ww, lsmask, nymin, mm, ndlon, nmonth, leap);
      }
      // do calculation on band (either over half band or just over one column dependeing on edgeflag)
      autocorr(lsmask, &ww, distances, lh[nmonth], mm, nxmin, nxmax, nymin, leap, ndlat, ndlon, nmonth, edgeflag);
    }
    dfree3(ww.data, (size_t) nlat, (size_t) nlon);
  }

  // interpolate lh and save to netcdf file
  savelh(lh, time, nxmin, nxmax, nymin, nymax, nmonth, leap);
  dfree3(lh, 14, (size_t) nlat);

  if (DBGFLG>2) {printf("  wavelength: return to main\n");fflush(NULL);}

  // free arrays
  dfree3(distances, (size_t) (nymax-nymin), (size_t) iterbounds[0]);
}

