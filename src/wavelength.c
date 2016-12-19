//
// Created by georg on 14/12/16.
//

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <math.h>
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
  sprintf(filepath, OUTPATH, nmonth);

  // count indicees
  start[0] = 0;
  count[0] = days[nmonth] * 24;

  // open file
  if ((retval = nc_open(filepath, NC_NOWRITE, &ncID))) ERR(retval);

  // get varID of vertical velocity
  if ((retval = nc_inq_varid(ncID, ZVEL, &varID))) ERR(retval);

  // read slices
  for (nn=nlatmin; nn<nlatmin + ww->nlat; nn++){
    for (mm=nlon; mm < nlon + 2 * ndlon; mm++) {
      // read slice
      if (mm < 0) {
        start[2] = (size_t) (lsmask->nlon + mm);
      } else if (mm > lsmask->nlon) {
        start[2] = (size_t) (mm - lsmask->nlon);
      }
      if (lsmask->data[nn][mm]==0.0){
        if ((retval = nc_get_vara_double(ncID, varID, start, count, &(ww->data)[0][nn-nlatmin][mm - nlon])))
        ERR(retval);
      }
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
  size_t start[3] = {0, 0, 0}, count[3]={days[nmonth] * 24, 1, 1}; // time, lats, lons
  char filepath[MAXCHARLEN];

  // shift data by one
  for (nn=0; nn < ww->nlat; nn++){
    for (mm=0; mm < 2 * ndlon - 1; mm++){
      for (nt=0; nt < days[nmonth]*24; nt++){
        ww->data[nt][nn][mm] = ww->data[nt][nn][mm+1];
      }
    }
  }

  // load one slice of data
  // set file to load from
  sprintf(filepath, OUTPATH, nmonth);

  // open file
  if ((retval = nc_open(filepath, NC_NOWRITE, &ncID))) ERR(retval);

  // get varID of vertical velocity
  if ((retval = nc_inq_varid(ncID, ZVEL, &varID))) ERR(retval);

  // read arrays
  for (nn=nlatmin; nn<nlatmin + ww->nlat; nn++){
    mm = nlon + 2 * ndlon;
    // read slice
    if (mm < 0) {
      start[2] = (size_t) (lsmask->nlon + mm);
    } else if (mm > lsmask -> nlon) {
      start[2] = (size_t) (mm - lsmask->nlon);
    }
    if (lsmask->data[nn][mm]==0.0){
      if ((retval = nc_get_vara_double(ncID, varID, start, count, &(ww->data)[0][nn-nlatmin][mm - nlon])))
      ERR(retval);
    }

  }

  // close nc file
  if ((retval = nc_close(ncID))) ERR(retval);
}

void autocorr(dat2d *lsmask, dat3d *ww, double ***distances, double **lh,
              int nlon, int nlonmin, int nlonmax, int nlatmin,
              int leap, int ndlat, int nmonth, int edgeflag) {
  // declare used variables
  size_t days[12] = {31, 28 + (size_t) leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int nn, nx, ny, nt, index, iterbounds[2];
  double *distance, *corr, *norm, avg;

  distance = dalloc(distance, CORRLEN);
  corr = dalloc(corr, CORRLEN);
  norm = dalloc(norm, CORRLEN);

//  for (index=0; index<CORRLEN; index++) distance[index]

  if (edgeflag==1) {
    // do calculation at the side determied by position nlon
  } else {
    for (nn=nlatmin; nn < nlatmin + ww->nlat; nn++) {
      // initialize wavelength with zero
      lh[nn-nlatmin][nlon-nlonmin] = 0.0;
      // boundaries for iteration over latitudes
      iterbounds[0] = (int) fmax(nn - ndlat/2, nlatmin);
      iterbounds[1] = (int) fmin(nn + ndlat/2, nlatmin + ww->nlat);
      // iterate over all time steps
      for(nt=0; nt<days[nmonth]; nt++) {
        // get average vertical velocity
        avg = davg2(ww->data[nt], ww->nlat, ww->nlon);
        // set correlation array and norms to zero
        for (index=0; index<CORRLEN; index++) corr[index] = norm[index] = 0.0;
        // iterate over lats and lons and add element to corresponding element in correlation array
        // note that we neglect any factor constant at the point since we evaluate maxima only
        for (ny=iterbounds[0]; ny<iterbounds[1]; ny++) {
          for (nx = 0; nx < ww->nlon; nx++) {
            // check distance and get index
            index = (int) ((distances[nn][ny][nx] - CORRMIN) / (CORRMAX - CORRMIN)*CORRLEN + 0.5);
            // add deviation to corr array and increase norm by 1
            corr[index] += (ww->data[nt][ny - nlatmin][nx] - avg) * (ww->data[nt][nn-nlatmin][nlon-nlonmin] - avg);
            norm[index]++;
          }
        }
        // normalize correlation with number of points in bin
        for (index=0; index<CORRLEN; index++) corr[index] /= norm[index];
        // get wavelegth from correlation array of point and add result to lh[nn-nlatmin][nlon-nlonmin]
        lh[nn-nlatmin][nlon-nlonmin] += dxmax(distance, corr, CORRLEN);
      }
      // normalize actual point in lh
      lh[nn-nlatmin][nlon-nlonmin] /= days[nmonth] * 24;
    }
  }

  free(distance);
  free(corr);
  free(norm);
}

void wavelength(dat2d *lsmask, int nxmin, int nxmax, int nymin, int nymax, int leap){
  size_t days[12] = {31, 28 + (size_t) leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  dat3d ww;
  size_t ntime;
  int nmonth, nlat, nlon, ndlon, nn, nx, ny, mm, itermax, edgeflag, ndlat;
  double latmax, dx;
  double ***distances, **lh;

  // prepare constants
  nlat = nymax - nymin; // number of points in latuitude of w-slice
  latmax = fmax(fabs(lsmask->lat[nymin]), fabs(lsmask->lat[nymax])); // maximum latitude
  dx = (lsmask->lon[1] - lsmask->lon[0])*EARTHRADIUS*cos(DEG2RAD(latmax)); // highest x-resolution in domain
  ndlon = (int) (CORRMAX/dx) + 2; // minimum number of points needed to include length from CORRMAX
  ndlat = 2 * (int) ((CORRMAX / fabs(lsmask->lat[1] - lsmask->lat[0])) + 2);

  // if domain is smaller than needed give error and abort calculation
  if ( 2 * ndlon > nxmax - nxmin) VALERR(2);
  nlon = 2 * ndlon; // number of points in lon of w-slice

  // check for repeating domain ('around the globe')
  edgeflag = 1;
  itermax = nxmax - nlon;
  if ((nxmax - nxmin) == lsmask->nlon) {
    itermax = lsmask->nlon;
    edgeflag = 0;
  };

  // set up distances arrays
  distances = dalloc3(distances, (size_t) (nymax-nymin), (size_t) ndlat, (size_t) ndlon);
  for (nn=0; nn<nymax-nymin; nn++){
    for (ny=0; ny<ndlat; ny++){
      for (nx=0; nx<ndlon; nx++){
        distances[nn][ny][nx] = dist(lsmask->lon[nx + nxmin], lsmask->lon[ndlon/2 + nxmin],
                                     lsmask->lat[ny + nymin], lsmask->lat[ny + nn]);
      }
    }
  }

  // iterate over months
  for (nmonth=0; nmonth<12; nmonth++) {
    ntime = days[nmonth] * 24;
    ww.data = dalloc3(ww.data, ntime, (size_t) nlat, (size_t) nlon);
    ww.nlat = nlat;
    ww.nlon = nlon;
    lh = dalloc2(lh, (size_t) nlat, (size_t) (nxmax-nxmin));
    // iterate over longitudes (+ handle edges of domains)
    for (mm=nxmin; mm < itermax; mm++){
      // load vertical velocity band
      if (mm==nxmin){
        loadw(&ww, lsmask, nymin, mm, ndlon, nmonth, leap);
      } else {
        advancew(&ww, lsmask, nymin, mm, ndlon, nmonth, leap);
      }
      // do calculation on band (either over half band or just over one column dependeing on edgeflag)
      autocorr(lsmask, &ww, distances, lh, mm, nxmin, nxmax, nymin, leap, ndlat, nmonth, edgeflag);
      // save data to netcdf file
    }

    dfree3(ww.data, ntime, (size_t) nlat);
    dfree2(lh, (size_t) nlat);
  }

  // free arrays
  dfree3(distances, (size_t) (nymax-nymin), (size_t) ndlat);
}

