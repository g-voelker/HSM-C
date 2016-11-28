//
// Created by georg on 24/11/16.
//
#include <stdlib.h>
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include <math.h>
#include "../lib/structs.h"
#include "../lib/dalloc.h"
#include "../lib/constants.h"
#include "header.h"
#include "input.h"

double distance(double lon1, double lon2, double lat1, double lat2){
  double dd;
  dd = EARTHRADIUS * 2 * acos(sin(DEG2RAD(lat1)) * sin(DEG2RAD(lat2)) +
                              cos(DEG2RAD(lat1)) * cos(DEG2RAD(lat2)) * cos(DEG2RAD(lon2 - lon1)));
  return(dd);
}

void divergence(dat2d *lsmask, int nxmin, int nxmax, int nymin, int nymax, int leap){
  // aux arrays
  double **uu, **vv, **ww, **mld, *aux, dys[nymax - nymin - 1], dx[nxmax - nxmin];

  // indexing
  int nn, mm, nmonth, nhour;
//  int days[12] = {31, 28 + leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  size_t days[12] = {31, 28 + (size_t) leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

  // netcdf
  int retval, ncID, varID[4], dimID[3];
  size_t start[3], count[3], latlen, lonlen;
  char filepath[MAXCHARLEN];
  size_t chunksize[3] = {28*24, CHUNK_LAT, CHUNK_LON};

  start[0] = start[1] = start[2] = 0;
  count[0] = 1;
  count[1] = (size_t) (nymax-nymin);
  count[2] = (size_t) (nxmax-nxmin);

  // allocate aux arrays
  uu = dalloc2(uu, (size_t) (nymax-nymin), (size_t) (nxmax-nxmin));
  vv = dalloc2(vv, (size_t) (nymax-nymin), (size_t) (nxmax-nxmin));
  ww = dalloc2(ww, (size_t) (nymax-nymin), (size_t) (nxmax-nxmin));
  mld = dalloc2(mld, (size_t) (nymax-nymin), (size_t) (nxmax-nxmin));
  aux = dalloc(aux, (size_t) ((nymax-nymin)*(nxmax-nxmin)));

  // calculate distance grid
  for (nn=nymin; nn<(nymax-1); nn++){ // on staggered grid
    dys[nn-nymin] = distance(0.0, 0.0, lsmask->lat[nn], lsmask->lat[nn+1]);
  }
  for (nn=nymin; nn<nymax; nn++){
    dx[nn-nymin] = distance(lsmask->lon[0], lsmask->lon[1], lsmask->lat[nn], lsmask->lat[nn]);
  }

  // loop over snap shots
  for (nmonth=0; nmonth< 12; nmonth++){
    chunksize[0] = days[nmonth]*24;
    sprintf(filepath, OUTPATH, nmonth+1);
    // open netcdf file and read u, v, mld
    if ((retval = nc_open(filepath, NC_WRITE, &ncID))) ERR(retval);

    // get variable IDs
    if ((retval = nc_inq_varid(ncID, XVEL, &varID[0]))) ERR(retval);
    if ((retval = nc_inq_varid(ncID, YVEL, &varID[1]))) ERR(retval);
    if ((retval = nc_inq_varid(ncID, MLD, &varID[2]))) ERR(retval);

    // get dimension IDs
    if ((retval = nc_inq_dimid(ncID, TIME, &dimID[0]))) ERR(retval);
    if ((retval = nc_inq_dimid(ncID, LATS, &dimID[1]))) ERR(retval);
    if ((retval = nc_inq_dimid(ncID, LONS, &dimID[2]))) ERR(retval);

    // check if there is at least 3x3 points
    if(nmonth==0){
      if((retval = nc_inq_dimlen(ncID, dimID[1], &latlen))) ERR(retval);
      if((retval = nc_inq_dimlen(ncID, dimID[2], &lonlen))) ERR(retval);
      if ((latlen<3)&(lonlen<3)){
        DIMERR(2);
      }
    }

    // define vertical velocity in data file
    if ((retval = nc_def_var(ncID, ZVEL, NC_DOUBLE, 3, dimID, &varID[3]))) ERR(retval);
    if ((retval = nc_def_var_chunking(ncID, varID[3], NC_CHUNKED, chunksize))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[3], UNITS, strlen(MPS), MPS))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[3], LONGNAME, strlen(ZVEL_LONG), ZVEL_LONG))) ERR(retval);

    for (nhour=0; nhour< (days[nmonth] * 24); nhour++) {
      // read variables u, v and mld
      start[0] = (size_t) nhour;
      if ((retval = nc_get_vara_double(ncID, varID[0], start, count, &aux[0]))) ERR(retval);
      for(nn=0; nn<(nymax-nymin); nn++){
        for (mm=0; mm<(nxmax-nxmin); mm++){
          uu[nn][mm] = aux[(nn)*(nxmax-nxmin) + (mm)];
        }
      }
      if ((retval = nc_get_vara_double(ncID, varID[1], start, count, &aux[0]))) ERR(retval);
      for(nn=0; nn<(nymax-nymin); nn++){
        for (mm=0; mm<(nxmax-nxmin); mm++){
          vv[nn][mm] = aux[(nn)*(nxmax-nxmin) + (mm)];
        }
      }
      if ((retval = nc_get_vara_double(ncID, varID[2], start, count, &aux[0]))) ERR(retval);
      for(nn=0; nn<(nymax-nymin); nn++){
        for (mm=0; mm<(nxmax-nxmin); mm++){
          mld[nn][mm] = aux[(nn)*(nxmax-nxmin) + (mm)];
        }
      }

      // go over all points and calculate divergence
      for (nn=nymin+1; nn<nymax-1; nn++){
        for (mm=nxmin+1; mm<nxmax-1; mm++){
          // corresponding index: (nn-nymin)*(nxmax-nxmin) + (mm-nxmin)
          // set condition for points on boundaries -> require all neighbors to exist
          if (((lsmask->data[nn-1][mm])==0)&((lsmask->data[nn+1][mm])==0)&
              ((lsmask->data[nn][mm-1])==0)&((lsmask->data[nn][mm+1])==0)){
            // all clear
            ww[nn-nymin][mm-nxmin] = (0.5 * (uu[nn+1-nymin][mm-nxmin] - uu[nn-nymin][mm-nxmin])/dys[nn-nymin] +
                    0.5 * (uu[nn-nymin][mm-nxmin] - uu[nn-1-nymin][mm-nxmin])/dys[nn-nymin-1] +
                    0.5 * (vv[nn-nymin][mm+1-nxmin] - vv[nn-nymin][mm-nxmin])/dx[nn-nymin] +
                    0.5 * (vv[nn-nymin][mm-nxmin] - vv[nn-nymin][mm-1-nxmin])/dx[nn-nymin]) * mld[nn-nymin][mm-nxmin];
          } else {
            ww[nn-nymin][mm-nxmin] = NC_FILL_DOUBLE;
          }
        }
      }
      for (nn=0; nn<(nymax-nymin); nn++){
        ww[nn][0] = ww[nn][nxmax-nxmin-1] = NC_FILL_DOUBLE;
      }
      for (mm=0; mm<(nxmax-nxmin); mm++){
        ww[0][mm] = ww[nymax-nymin-1][mm] = NC_FILL_DOUBLE;
      }
//      if ((nmonth==0)&(nhour==0)) {
//        for (nn=0; nn<3; nn++){
//          printf("%e\t%e\t%e\n", ww[nn][0], ww[nn][1], ww[nn][2]);
//          fflush(NULL);
//        }
//      }
      // dump w to file
      for(nn=0; nn<(nymax-nymin); nn++) {
        for (mm = 0; mm < (nxmax - nxmin); mm++) {
          aux[(nn) * (nxmax - nxmin) + (mm)] = ww[nn][mm];
        }
      }
      if ((retval = nc_put_vara_double(ncID, varID[3], start, count, &aux[0]))) ERR(retval);
    }
    // close netcdf file
    if ((retval = nc_close(ncID))) ERR(retval);
  }

  // free temporary arrays
  dfree2(uu, (size_t) (nymax - nymin));
  dfree2(vv, (size_t) (nymax - nymin));
  dfree2(ww, (size_t) (nymax - nymin));
  dfree2(mld, (size_t) (nymax - nymin));
  free(aux);
}