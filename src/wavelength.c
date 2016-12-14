//
// Created by georg on 14/12/16.
//

#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include "../lib/structs.h"
#include "input.h"
#include "header.h"
#include "save.h"


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
    for (mm=nlon-ndlon; mm < nlon+ndlon; mm++) {
      // read slice
      if (mm < 0) {
        start[2] = (size_t) (lsmask->nlon + mm);
      } else if (mm > lsmask -> nlon) {
        start[2] = (size_t) (mm - lsmask->nlon);
      }
      if (lsmask->data[nn][mm]==0.0){
        if ((retval = nc_get_vara_double(ncID, varID, start, count, &(ww->data)[0][nn-nlatmin][mm - (nlon - ndlon)])))
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
    for (mm=0; mm< 2*ndlon - 1; mm++){
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
    mm = nlon + ndlon;
    // read slice
    if (mm < 0) {
      start[2] = (size_t) (lsmask->nlon + mm);
    } else if (mm > lsmask -> nlon) {
      start[2] = (size_t) (mm - lsmask->nlon);
    }
    if (lsmask->data[nn][mm]==0.0){
      if ((retval = nc_get_vara_double(ncID, varID, start, count, &(ww->data)[0][nn-nlatmin][mm - (nlon - ndlon)])))
      ERR(retval);
    }

  }

  // close nc file
  if ((retval = nc_close(ncID))) ERR(retval);
}

void wavelength(dat2d *lsmask, int nxmin, int nxmax, int nymin, int nymax, int leap){
  // load longitude band and iterate over globe
}

void getlh(dat2d *lsmask, dat1d *lh, int nn, int mm){
  // get indicees of points to load

  // load matrix

  // analyse time step by time step

  // average time bins together

  // interpolate reduced time series

  // set values in struct

}