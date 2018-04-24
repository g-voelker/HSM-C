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
#include "../lib/macros.h"
#include "header.h"

double dist(double lon1, double lon2, double lat1, double lat2){
  double dd = 0.0;
  if ((lat1!=lat2)|(lon1!=lon2)){
    dd = EARTHRADIUS * 2 * acos(sin(DEG2RAD(lat1)) * sin(DEG2RAD(lat2)) +
                                cos(DEG2RAD(lat1)) * cos(DEG2RAD(lat2)) * cos(DEG2RAD(lon2 - lon1)));
  }
  return(dd);
}

void divergence(dat2d *lsmask, double *params, char **paths,
                int nxmin, int nxmax, int nymin, int nymax, int leap, int hemflag){
  // aux arrays
  double **uu, **vv, *aux, dys[nymax - nymin - 1], dx[nymax - nymin];
  double *vc, *vn, *vs, *uc, *ue, *uw, *ww, *mld;

  // indexing
  int nn, mm, nmonth, nhour;
  size_t days[14] = {31, 31, 28 + (size_t) leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31};

  // netcdf
  int retval, ncID = 0, varID[4], dimID[3];
  size_t start[3], count[3], latlen, lonlen;
  char filepath[MAXCHARLEN];
  size_t chunksize[3] = {28*24, CHUNK_LAT, CHUNK_LON};

  start[0] = start[1] = start[2] = 0;
  count[0] = count[1] = count[2] = 1;

  if (DBGFLG>2) {printf("  divergence: calculate distances in grid\n");fflush(NULL);}
  // calculate distance grid
  for (nn=nymin; nn<(nymax-1); nn++){ // on staggered grid
    printf("nymin: %d", nymin);fflush(NULL);
    printf("nymax: %d", nymax);fflush(NULL);
    printf("nn: %d", nn);fflush(NULL);
    printf("lat[nn]: %f", lsmask->lat[nn]);fflush(NULL);
    printf("lat[nn]: %f", lsmask->lat[nn+1]);fflush(NULL);
    dys[nn-nymin] = dsgn(lsmask->lat[nn+1] - lsmask->lat[nn]) * dist(0.0, 0.0, lsmask->lat[nn], lsmask->lat[nn+1]);
  }
  for (nn=nymin; nn<nymax; nn++){
    dx[nn-nymin] = dsgn(lsmask->lon[1] - lsmask->lon[0]) * dist(lsmask->lon[0], lsmask->lon[1], lsmask->lat[nn], lsmask->lat[nn]);
  }

  if (DBGFLG>2) {printf("  divergence: loop over all months and snapshots\n");fflush(NULL);}
  // loop over snapshots
  for (nmonth=0; nmonth< 14; nmonth++){
    // set count and chunksize in time
    count[0] = chunksize[0] = days[nmonth]*24;


    // allocate auxiliary arrays
    vc = dalloc(vc, days[nmonth]*24);
    vn = dalloc(vn, days[nmonth]*24);
    vs = dalloc(vs, days[nmonth]*24);
    uc = dalloc(uc, days[nmonth]*24);
    ue = dalloc(ue, days[nmonth]*24);
    uw = dalloc(uw, days[nmonth]*24);
    ww = dalloc(ww, days[nmonth]*24);
    mld = dalloc(mld, days[nmonth]*24);

    // note: the order of computation doesn't matter here
    if (hemflag==0) { // northern hemisphere
      if (nmonth == 0) {
        sprintf(filepath, paths[3], (int) params[10],  1);
      } else if (nmonth == 13) {
        sprintf(filepath, paths[3], (int) params[10], 12);
      } else {
        sprintf(filepath, paths[5], (int) params[10], nmonth);
      }
    } else if (hemflag==1) { // southern hemisphere
      if (nmonth == 0) {
        sprintf(filepath, paths[4], (int) params[10], 1);
      } else if (nmonth == 13) {
        sprintf(filepath, paths[4], (int) params[10], 12);
      } else {
        sprintf(filepath, paths[4], (int) params[10], nmonth);
      }
    } else { // if on no known hemisphere throw error
      GENERR
    }

    // open netcdf file and read u, v, mld
    if ((retval = nc_open(filepath, NC_WRITE, &ncID))) ERR(retval);

    // get variable IDs
    if ((retval = nc_inq_varid(ncID, XVEL, &varID[0]))) ERR(retval);
    if ((retval = nc_inq_varid(ncID, YVEL, &varID[1]))) ERR(retval);
    if ((retval = nc_inq_varid(ncID, MLD, &varID[2]))) ERR(retval);
    if ((retval = nc_inq_varid(ncID, ZVEL, &varID[3]))) ERR(retval);

    // get dimension IDs
    if ((retval = nc_inq_dimid(ncID, TIME, &dimID[0]))) ERR(retval);
    if ((retval = nc_inq_dimid(ncID, LATS, &dimID[1]))) ERR(retval);
    if ((retval = nc_inq_dimid(ncID, LONS, &dimID[2]))) ERR(retval);

    // check if there is at least 3x3 points
    if(nmonth==0){
      if((retval = nc_inq_dimlen(ncID, dimID[1], &latlen))) ERR(retval);
      if((retval = nc_inq_dimlen(ncID, dimID[2], &lonlen))) ERR(retval);
      if ((latlen<3) | (lonlen<3)){
        DIMERR;
      }
    }

    // iterate over lats and lons
    printf("nmonth: %d\n", nmonth);fflush(NULL);
    for (nn = nymin+1; nn < nymax-1; nn++){
//      printf("nlat: %d\n", nn);fflush(NULL);
      for (mm = nxmin+1; mm < nxmax-1; mm++){
//        if (nn==271) printf("nlon: %d\n", mm);fflush(NULL);
        if (lsmask->data[nn-1][mm] +
            lsmask->data[nn+1][mm] +
            lsmask->data[nn][mm-1] +
            lsmask->data[nn][mm+1] +
            lsmask->data[nn][mm] == 0.0) {
          // read adjacent and central time series from netcdf file
          start[1] = (size_t) nn - nymin;
          start[2] = (size_t) mm - nxmin - 1;
          if ((retval = nc_get_vara_double(ncID, varID[0], start, count, &uw[0]))) ERR(retval);
          start[1] = (size_t) nn - nymin;
          start[2] = (size_t) mm - nxmin + 1;
          if ((retval = nc_get_vara_double(ncID, varID[0], start, count, &ue[0]))) ERR(retval);
          start[1] = (size_t) nn - nymin + 1;
          start[2] = (size_t) mm - nxmin;
          if ((retval = nc_get_vara_double(ncID, varID[1], start, count, &vn[0]))) ERR(retval);
          start[1] = (size_t) nn - nymin - 1;
          start[2] = (size_t) mm - nxmin;
          if ((retval = nc_get_vara_double(ncID, varID[1], start, count, &vs[0]))) ERR(retval);
          start[1] = (size_t) nn - nymin;
          start[2] = (size_t) mm - nxmin;
          if ((retval = nc_get_vara_double(ncID, varID[0], start, count, &uc[0]))) ERR(retval);
          if ((retval = nc_get_vara_double(ncID, varID[1], start, count, &vc[0]))) ERR(retval);
          if ((retval = nc_get_vara_double(ncID, varID[2], start, count, &mld[0]))) ERR(retval);

          // calculate vertical velocity
          for (nhour = 0; nhour < count[0]; nhour++) {
            ww[nhour] = 0.5 * ((vn[nhour] - vc[nhour]) / dys[nn - nymin] +
                               (vc[nhour] - vs[nhour]) / dys[nn - nymin - 1] +
                               (ue[nhour] - uc[nhour]) / dx[nn - nymin - 1] +
                               (uc[nhour] - uw[nhour]) / dx[nn - nymin - 1]) * mld[nhour];
          }
        } else {
          for (nhour = 0; nhour < count[0]; nhour++) {
            ww[nhour] = NC_FILL_DOUBLE;
          }
        }

        // dump point to netCDF file
        if ((retval = nc_put_vara_double(ncID, varID[3], start, count, &ww[0]))) ERR(retval);
      }
    }
    // close file
    if ((retval = nc_close(ncID))) ERR(retval);

    // free arrays for next loop
    free(vn);
    free(vs);
    free(vc);
    free(ue);
    free(uw);
    free(uc);
    free(ww);
    free(mld);
  }
  if (DBGFLG>2) {printf("  divergence: return to main\n");fflush(NULL);}
}