//
// Created by georg on 21/11/16.
//

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <netcdf.h>
#include "../lib/ialloc.h"
#include "../lib/structs.h"
#include "../lib/slicer.h"
#include "save.h"
#include "input.h"
#include "header.h"
#include "../lib/dalloc.h"
#include "../lib/macros.h"

int sum(size_t *array, int num){
  int nn;
  size_t sum=0;
  for (nn=0; nn<num; nn++){
    sum += array[nn];
  }
  return((int) sum);
}

void initnc_write(dat2d *lsmask, int* time, int nxmin, int nxmax, int nymin, int nymax, int leap, int splitflag){
  int *timeSliceS, *timeSliceH, varID[7], nn = 0;
  size_t nt;
  size_t chunksize[3] = {28*24, CHUNK_LAT, CHUNK_LON};
  int ncID = 0, retval, dimID[3], dimVarID[3];
  char filepath[MAXCHARLEN];
  size_t days[14] = {31, 31, 28 + (size_t) leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31};
  double *latSlice, *lonSlice;

  if (DBGFLG>2) {printf("  initnc: slice lats and lons to interval\n"); fflush(NULL);}

  latSlice = dalloc(latSlice, (size_t) (nymax - nymin));
  latSlice = dslice1(lsmask->lat, latSlice, nymin, nymax);

  lonSlice = dalloc(lonSlice, (size_t) (nxmax - nxmin));
  lonSlice = dslice1(lsmask->lon, lonSlice, nxmin, nxmax);

  if (DBGFLG>2) {printf("  initnc: init nc files\n"); fflush(NULL);}
  // set reference time for time axis

  for (nn=0; nn < 14; nn++) {
    chunksize[0] = days[nn]*24;
    // set filename
    if (splitflag==0) {
      if (nn == 0) {
        sprintf(filepath, AUXPATH_N, 1);
      } else if (nn == 13) {
        sprintf(filepath, AUXPATH_N, 12);
      } else {
        sprintf(filepath, OUTPATH_N, nn);
      }
    } else if (splitflag==1) {
      if (nn == 0) {
        sprintf(filepath, AUXPATH_S, 1);
      } else if (nn == 13) {
        sprintf(filepath, AUXPATH_S, 12);
      } else {
        sprintf(filepath, OUTPATH_S, nn);
      }
    }
    // created the file with unlimited dimension
    if ((retval = nc_create(filepath, NC_NETCDF4 | NC_CLOBBER, &ncID))) ERR(retval);

    // set dimensions
    // time
    if ((retval = nc_def_dim(ncID, TIME, days[nn] * 24, &dimID[0]))) ERR(retval);
    if ((retval = nc_def_var(ncID, TIME, NC_INT, 1, &dimID[0], &dimVarID[0]))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, dimVarID[0], UNITS, strlen(HOURS), HOURS))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, dimVarID[0], LONGNAME, strlen(TIME), TIME))) ERR(retval);

    // latitude
    if ((retval = nc_def_dim(ncID, LATS, (size_t) (nymax - nymin), &dimID[1]))) ERR(retval);
    if ((retval = nc_def_var(ncID, LATS, NC_DOUBLE, 1, &dimID[1], &dimVarID[1]))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, dimVarID[1], UNITS, strlen(DEGREES_NORTH), DEGREES_NORTH))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, dimVarID[1], LONGNAME, strlen(LATS), LATS))) ERR(retval);
    // longitude
    if ((retval = nc_def_dim(ncID, LONS, (size_t) (nxmax - nxmin), &dimID[2]))) ERR(retval);
    if ((retval = nc_def_var(ncID, LONS, NC_DOUBLE, 1, &dimID[2], &dimVarID[2]))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, dimVarID[2], UNITS, strlen(DEGREES_EAST), DEGREES_EAST))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, dimVarID[2], LONGNAME, strlen(LONS), LONS))) ERR(retval);

    // define matrices with Fill Values that are replaced later
    // u component
    if ((retval = nc_def_var(ncID, XVEL, NC_DOUBLE, 3, dimID, &varID[0]))) ERR(retval);
    if ((retval = nc_def_var_chunking(ncID, varID[0], NC_CHUNKED, chunksize))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[0], UNITS, strlen(MPS), MPS))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[0], LONGNAME, strlen(XVEL_LONG), XVEL_LONG))) ERR(retval);
    // v component
    if ((retval = nc_def_var(ncID, YVEL, NC_DOUBLE, 3, dimID, &varID[1]))) ERR(retval);
    if ((retval = nc_def_var_chunking(ncID, varID[1], NC_CHUNKED, chunksize))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[1], UNITS, strlen(MPS), MPS))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[1], LONGNAME, strlen(YVEL_LONG), YVEL_LONG))) ERR(retval);
    // define vertical velocity in data file
    if ((retval = nc_def_var(ncID, ZVEL, NC_DOUBLE, 3, dimID, &varID[2]))) ERR(retval);
    if ((retval = nc_def_var_chunking(ncID, varID[2], NC_CHUNKED, chunksize))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[2], UNITS, strlen(MPS), MPS))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[2], LONGNAME, strlen(ZVEL_LONG), ZVEL_LONG))) ERR(retval);
    // mixed layer depth
    if ((retval = nc_def_var(ncID, MLD, NC_DOUBLE, 3, dimID, &varID[3]))) ERR(retval);
    if ((retval = nc_def_var_chunking(ncID, varID[3], NC_CHUNKED, chunksize))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[3], UNITS, strlen(METER), METER))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[3], LONGNAME, strlen(MLD_LONG), MLD_LONG))) ERR(retval);
    // wind work
    if ((retval = nc_def_var(ncID, EIN, NC_DOUBLE, 3, dimID, &varID[4]))) ERR(retval);
    if ((retval = nc_def_var_chunking(ncID, varID[4], NC_CHUNKED, chunksize))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[4], UNITS, strlen(WPM2), WPM2))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[4], LONGNAME, strlen(EIN_LONG), EIN_LONG))) ERR(retval);
    // energy flux radiated
    if ((retval = nc_def_var(ncID, EOUT, NC_DOUBLE, 3, dimID, &varID[5]))) ERR(retval);
    if ((retval = nc_def_var_chunking(ncID, varID[5], NC_CHUNKED, chunksize))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[5], UNITS, strlen(WPM2), WPM2))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[5], LONGNAME, strlen(EOUT_LONG), EOUT_LONG))) ERR(retval);
    // horizontal length scale
    if ((retval = nc_def_var(ncID, LH, NC_DOUBLE, 3, dimID, &varID[6]))) ERR(retval);
    if ((retval = nc_def_var_chunking(ncID, varID[5], NC_CHUNKED, chunksize))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[5], UNITS, strlen(METER), METER))) ERR(retval);
    if ((retval = nc_put_att_text(ncID, varID[5], LONGNAME, strlen(LH_LONG), LH_LONG))) ERR(retval);

    // this is the end of the definition mode

    // allocate time array with variable length and slice from time axis
    timeSliceS = ialloc(timeSliceS, days[nn] * 24);
    timeSliceH = ialloc(timeSliceH, days[nn] * 24);
    timeSliceS = islice1(time, timeSliceS, sum(days, nn) * 24, sum(days, nn+1) * 24);
    for (nt=0; nt<days[nn]*24; nt++){
      timeSliceH[nt] = timeSliceS[nt] / 3600 + 613608;
    }
    // save time slice
    if ((retval = nc_put_var_int(ncID, dimVarID[0], &timeSliceH[0]))) ERR(retval);
    // free time array for next time month
    free(timeSliceS);
    free(timeSliceH);

    // save latitude
    if ((retval = nc_put_var_double(ncID, dimVarID[1], &latSlice[0]))) ERR(retval);

    // save longitude
    if ((retval = nc_put_var_double(ncID, dimVarID[2], &lonSlice[0]))) ERR(retval);

    // close netcdf file
    if ((retval = nc_close(ncID))) ERR(retval);
  }

  // free sliced arrays
  free(latSlice);
  free(lonSlice);
}

void initnc(dat2d *lsmask, int* time, int nxmin, int nxmax, int nymin, int nymax, int leap, int *nlat5, int *slat5){
  int splitflag = 0, nn;

  *nlat5 = NC_FILL_INT;
  *slat5 = NC_FILL_INT;

  // check if we are operating on northern / southern / both hemispheres
  if ((lsmask->lat[nymin]>0) & (lsmask->lat[nymax]>0)) {
    splitflag = 0;
  } else if ((lsmask->lat[nymin]<0)&(lsmask->lat[nymax]<0)) {
    splitflag = 1;
  } else {
    splitflag = 2;
    // throw error when region is invalid (within 5 degr. from equator
    if ((dabs(lsmask->lat[nymin]) < 5) & (dabs(lsmask->lat[nymax]) < 5)){
      DOMAINERR;
    }
  }

  // write accordingly
  if (splitflag==0){
    initnc_write(lsmask, time, nxmin, nxmax, nymin, nymax, leap, splitflag);
  } else if (splitflag==1) {
    initnc_write(lsmask, time, nxmin, nxmax, nymin, nymax, leap, splitflag);
  } else if (splitflag==2) {
    // find nlat5 and slat5
    for (nn=nymin; nn<nymax; nn++){
      if (dabs(lsmask->lat[nn]) < 5) break;
    }

    if (lsmask->lat[nn] > 0){
      memcpy(nlat5, &nn, sizeof(int));
    } else {
      memcpy(slat5, &nn, sizeof(int));
    }

    for (nn=nymax; nn>nymin; nn--){
      if (dabs(lsmask->lat[nn]) < 5) break;
    }

    if (lsmask->lat[nn] > 0){
      memcpy(nlat5, &nn, sizeof(int));
    } else {
      memcpy(slat5, &nn, sizeof(int));
    }

    // if nlat5 or slat5 did not receive a value throw an error
    if ((*nlat5==NC_FILL_INT)|(*slat5==NC_FILL_INT)){
      GENERR
    }

    // write two files with according nymin / nymax and according splitflag
    if (nymin == *nlat5) {
      // in this case the border at nymin is within 5 degr of the  the equator (NORTH)
      initnc_write(lsmask, time, nxmin, nxmax, *slat5, nymax, leap, 1);
    } else if (nymax == *nlat5) {
      // in this case the border at nymax is within 5 degr of the  the equator (NORTH)
      initnc_write(lsmask, time, nxmin, nxmax, nymin, *slat5, leap, 1);
    } else if (nymin == *slat5) {
      // in this case the border at nymin is within 5 degr of the  the equator (SOUTH)
      initnc_write(lsmask, time, nxmin, nxmax, *nlat5, nymax, leap, 0);
    } else if (nymax == *slat5) {
      // in this case the border at nymax is within 5 degr of the  the equator (SOUTH)
      initnc_write(lsmask, time, nxmin, nxmax, nymin, *nlat5, leap, 0);
    } else {
      // in this case there is a full domain that stretches on both sides of the equator
      if (lsmask->lat[nymin]<0) {

        initnc_write(lsmask, time, nxmin, nxmax, nymin, *slat5, leap, 1);
        initnc_write(lsmask, time, nxmin, nxmax, *nlat5, nymax, leap, 0);

      } else if (lsmask->lat[nymin]>0) {

        initnc_write(lsmask, time, nxmin, nxmax, *slat5, nymax, leap, 1);
        initnc_write(lsmask, time, nxmin, nxmax, nymin, *nlat5, leap, 0);

      }
    }
  }
}

void savePoint(dat2d *lsmask, double* uu, double* vv, double* mld, double* taux, double* tauy,
               int* time, int nxmin, int nx, int nymin, int nymax, int ny, int leap, int nlat5, int slat5){
  size_t days[14] = {31, 31, 28 + (size_t) leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31};
  int nn = 0, nt = 0, mon = 0, ncID, retval, varID, mm = 0;
  int *index;
  int iterbounds[2];
  size_t start[3], count[3];
  time_t aux;
  struct tm *auxTime;
  double *uSlice, *vSlice, *mldSlice, *windWork;
  char filepath[MAXCHARLEN];

  for (nn=0; nn<14; nn++) {
    if (DBGFLG>2) {printf("  savePoint: processing month %d\n", nn); fflush(NULL);}
    if (DBGFLG>2) {printf("  savePoint: slicing data\n"); fflush(NULL);}
    // set slicing indices according to month
    iterbounds[0] = sum(days, nn) * 24;
    iterbounds[1] = sum(days, nn+1) * 24;

    // slice data according to month
    uSlice = dalloc(uSlice, (size_t) (iterbounds[1] - iterbounds[0]));
    vSlice = dalloc(vSlice, (size_t) (iterbounds[1] - iterbounds[0]));
    mldSlice = dalloc(mldSlice, (size_t) (iterbounds[1] - iterbounds[0]));
    windWork = dalloc(windWork, (size_t) (iterbounds[1] - iterbounds[0]));


    uSlice = dslice1(uu, uSlice, iterbounds[0], iterbounds[1]);
    vSlice = dslice1(vv, vSlice, iterbounds[0], iterbounds[1]);
    mldSlice = dslice1(mld, mldSlice, iterbounds[0], iterbounds[1]);

    if (DBGFLG>2) {printf("  savePoint: calculating wind work\n"); fflush(NULL);}
    for (mm=iterbounds[0]; mm < iterbounds[1]; mm++){
      windWork[mm - iterbounds[0]] = taux[mm] * uSlice[mm - iterbounds[0]] + tauy[mm] * vSlice[mm - iterbounds[0]];
    }

    // set filepath depending on hemisphere
    if (lsmask->lat[ny]>5) {
      if (nn == 0) {
        sprintf(filepath, AUXPATH_N, 1);
      } else if (nn == 13) {
        sprintf(filepath, AUXPATH_N, 12);
      } else {
        sprintf(filepath, OUTPATH_N, nn);
      }
    } else if (lsmask->lat[ny]<-5) {
      if (nn == 0) {
        sprintf(filepath, AUXPATH_S, 1);
      } else if (nn == 13) {
        sprintf(filepath, AUXPATH_S, 12);
      } else {
        sprintf(filepath, OUTPATH_S, nn);
      }
    }

    // set start index according to hemisphere
    if (((lsmask->lat[nymin]>0)&(lsmask->lat[nymax]>0))|
        ((lsmask->lat[nymin]<0)&(lsmask->lat[nymax]<0))|
        ((lsmask->lat[ny] > 5)&(lsmask->lat[nymin] > lsmask->lat[nymax]))|
        ((lsmask->lat[ny] > -5)&(lsmask->lat[nymin] < lsmask->lat[nymax]))){

      start[1] = (size_t) (ny-nymin);

    } else if ((lsmask->lat[ny] > 5)&(lsmask->lat[nymin] < lsmask->lat[nymax])) {

      start[1] = (size_t) (ny-nlat5);

    } else if ((lsmask->lat[ny] < -5)&(lsmask->lat[nymin] > lsmask->lat[nymax])) {

      start[1] = (size_t) (ny-slat5);

    } else {
      // if lat of point is in forbidden latitude band throw an error
      GENERR
    }

    // set hyperslab indicees
    start[0] = (size_t) 0;
    count[0] = (size_t) (iterbounds[1] - iterbounds[0]);
    count[1] =  1;
    start[2] = (size_t) (nx-nxmin);
    count[2] =  1;

//    printf("%d\n", ny);
//    printf("%d\n", slat5);
//    printf("%d, %d, %d, %d, %d, %d\n",
//           (int) start[0], (int) count[0],
//           (int) start[1], (int) count[1],
//           (int) start[2], (int) count[2]);

    if (DBGFLG>2) {printf("  savePoint: saving data to nc file\n"); fflush(NULL);}
    // open nc file
    if ((retval = nc_open(filepath, NC_WRITE, &ncID))) ERR(retval);

    // get variable ID and write hyperslab into file
    if ((retval = nc_inq_varid(ncID, XVEL, &varID))) ERR(retval);
    if ((retval = nc_put_vara_double(ncID, varID, start, count, &uSlice[0]))) ERR(retval);

    // get variable ID and write hyperslab into file
    if ((retval = nc_inq_varid(ncID, YVEL, &varID))) ERR(retval);
    if ((retval = nc_put_vara_double(ncID, varID, start, count, &vSlice[0]))) ERR(retval);

    // get variable ID and write hyperslab into file
    if ((retval = nc_inq_varid(ncID, MLD, &varID))) ERR(retval);
    if ((retval = nc_put_vara_double(ncID, varID, start, count, &mldSlice[0]))) ERR(retval);

    // get variable ID and write hyperslab into file
    if ((retval = nc_inq_varid(ncID, EIN, &varID))) ERR(retval);
    if ((retval = nc_put_vara_double(ncID, varID, start, count, &windWork[0]))) ERR(retval);

    // close nc file
    if ((retval = nc_close(ncID))) ERR(retval);

    free(uSlice);
    free(vSlice);
    free(mldSlice);
    free(windWork);
  }
  if (DBGFLG>2) {printf("  savePoint: return to main\n"); fflush(NULL);}
}

void savelh(dat2d *lsmask, double ***lh, int *time,
            int nxmin, int nxmax, int nymin, int nymax, int leap, int hemflag){
  size_t days[12] = {31, 28 + (size_t) leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int nx, ny, nt, nint, ncID, retval, varID, nmonth = 0;
  size_t start[3], count[3];
  int *lhTime;
  struct tm tcon;
  double *lhSlice;
  char filepath[MAXCHARLEN];

  // set lh time for interpolation
  lhTime = ialloc(lhTime, 14);
  for (nt=0; nt<14; nt++){
    if (nt==0){
      tcon.tm_year = YEAR - 1901;
      tcon.tm_mon = 11;
    } else if (nt==13) {
      tcon.tm_year = YEAR - 1899;
      tcon.tm_mon = 0;
    } else {
      tcon.tm_year = YEAR - 1900;
      tcon.tm_mon = nt-1;
    }
    tcon.tm_mday = 15;
    tcon.tm_hour = 1;
    tcon.tm_min = 0;
    tcon.tm_sec = 0;
    // tcon.tm_gmtoff = 0;
    tcon.tm_isdst = 0;
    lhTime[nt] =  (int) mktime(&tcon);
  }

  if (DBGFLG>2) {printf("  savelh: interpolate and save wavelength to data file\n"); fflush(NULL);}

  // iterate over months
  for (nmonth=0; nmonth<12; nmonth++) {
    // allocate data according to month
    lhSlice = dalloc(lhSlice, (size_t) days[nmonth]*24);

    // set filepath as above
    if (hemflag==0) {
      sprintf(filepath, OUTPATH_N, nmonth+1);
    } else if (hemflag==1) {
      sprintf(filepath, OUTPATH_S, nmonth+1);
    }

    // set hyperslab indicees
    start[0] = (size_t) 0;
    count[0] = (size_t) (days[nmonth]*24);

    // open nc file and get variable ID
    if ((retval = nc_open(filepath, NC_WRITE, &ncID))) ERR(retval);
    if ((retval = nc_inq_varid(ncID, LH, &varID))) ERR(retval);

    // start iteration here
    for (ny=nymin; ny<nymax; ny++) {
      for (nx=nxmin; nx<nxmax; nx++) {
        // interpolate in time
        for (nt=0; nt<days[nmonth]*24; nt++) {
          // check position relative to monthly values
          for (nint = 1; nint < 14; nint++) {
            if (lhTime[nint] > time[nt]) {
              break;
            }
          }
          // linearly interpolate
          lhSlice[nt] = (lh[nint - 1][ny - nymin][nx - nxmin] +
                        (lh[nint][ny - nymin][nx - nxmin] - lh[nint - 1][ny - nymin][nx - nxmin]) /
                        (lhTime[nint] - lhTime[nint - 1]) *
                        (time[nt] - lhTime[nint - 1]));
        }

        // set hyperlsab
        start[1] = (size_t) (ny - nymin);
        count[1] = 1;
        start[2] = (size_t) (nx - nxmin);
        count[2] = 1;

        // write hyperslab into file
        if ((retval = nc_put_vara_double(ncID, varID, start, count, &lhSlice[0]))) ERR(retval);
      }
    }

    // close nc file
    if ((retval = nc_close(ncID))) ERR(retval);

    free(lhSlice);
  }
  free(lhTime);
  if (DBGFLG>2) {printf("  savelh: done.\n"); fflush(NULL);}
}

void savePointHybrid(dat2d *lsmask, dat1d *Eout, int ny, int nymin, int nymax, int nx, int nxmin, int leap, int nlat5, int slat5){
  size_t days[14] = {31, 28 + (size_t) leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int nmonth = 0, ncID, retval, varID;
  int iterbounds[2];
  size_t start[3], count[3];
  double *EoutSlice;
  char filepath[MAXCHARLEN];

  for (nmonth=0; nmonth<12; nmonth++) {
    if (DBGFLG>2) {printf("  savePoint: processing month %d\n", nmonth + 1); fflush(NULL);}
    if (DBGFLG>2) {printf("  savePoint: slicing data\n"); fflush(NULL);}
    // set slicing indices according to month
    iterbounds[0] = sum(days, nmonth) * 24;
    iterbounds[1] = sum(days, nmonth+1) * 24;

    // slice data according to month
    EoutSlice = dalloc(EoutSlice, (size_t) (iterbounds[1] - iterbounds[0]));
    EoutSlice = dslice1(Eout->data, EoutSlice, iterbounds[0], iterbounds[1]);

    // set filepath as above
    if (lsmask->lat[ny]>5) {
      sprintf(filepath, OUTPATH_N, nmonth+1);
    } else if (lsmask->lat[ny]<-5) {
      sprintf(filepath, OUTPATH_S, nmonth+1);
    } else {
      // if lat of point is in forbidden latitude band throw an error
      GENERR
    }

    // set start index according to hemisphere
    if (((lsmask->lat[nymin]>0)&(lsmask->lat[nymax]>0))|
        ((lsmask->lat[nymin]<0)&(lsmask->lat[nymax]<0))|
        ((lsmask->lat[ny] > 5)&(lsmask->lat[nymin] > lsmask->lat[nymax]))|
        ((lsmask->lat[ny] > -5)&(lsmask->lat[nymin] < lsmask->lat[nymax]))){

      start[1] = (size_t) (ny-nymin);

    } else if ((lsmask->lat[ny] > 5)&(lsmask->lat[nymin] < lsmask->lat[nymax])) {

      start[1] = (size_t) (ny-nlat5);

    } else if ((lsmask->lat[ny] < -5)&(lsmask->lat[nymin] > lsmask->lat[nymax])) {

      start[1] = (size_t) (ny-slat5);

    } else {
      // if lat of point is in forbidden latitude band throw an error
      GENERR
    }

    // set hyperslab indicees
    start[0] = (size_t) 0;
    count[0] = (size_t) (iterbounds[1] - iterbounds[0]);
    count[1] =  1;
    start[2] = (size_t) (nx-nxmin);
    count[2] =  1;

    if (DBGFLG>2) {printf("  savePoint: saving data to nc file\n"); fflush(NULL);}
    // open nc file
    if ((retval = nc_open(filepath, NC_WRITE, &ncID))) ERR(retval);

    // get variable ID and write hyperslab into file
    if ((retval = nc_inq_varid(ncID, EOUT, &varID))) ERR(retval);
    if ((retval = nc_put_vara_double(ncID, varID, start, count, &EoutSlice[0]))) ERR(retval);

    // close nc file
    if ((retval = nc_close(ncID))) ERR(retval);

    free(EoutSlice);
  }
  if (DBGFLG>2) {printf("  savePoint: return to main\n"); fflush(NULL);}

}