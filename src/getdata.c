#include <netcdf.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../lib/dalloc.h"
#include "../lib/ialloc.h"
#include "../lib/structs.h"
#include "header.h"
#include "save.h"
// #include <grib_api.h>
// #include <unistd.h> (sleep)


void correct(double *input, int length){
  /* data arrays in grib files are averages over periodes
   * 0-1
   * 0-2
   * ...
   * 0-6
   * to correct to averages of the respective hour we correct with
   * input(n)*(n%6+1)-input(n-1)*(n%6)
   * 
   * inputs: double *array, int lengthOfARray
   * corrects the array in place.
   */
  
  double *temp;

  temp = dalloc(temp, (size_t) length);
  int nn;

  for (nn=0; nn<length; nn++){
    temp[nn] = input[nn];
  }

  for (nn=1; nn<length; nn++) {
    input[nn] = temp[nn]*(nn%6+1) - temp[nn-1]*(nn%6);
  }

  free(temp);
}

/* data to get, given nlat, nlon:
 * mld (in m)
 * tau_x (in N/m^2)
 * tau_y (in N/m^2)
 * time axis (in seconds from 1/1/1970 00:00:00)
 * 
 * modifications that are done:
 * interpolate mld to given year's timeline
 * correct stress values as they come in averages
 */

void getdata(double *params, char **paths,
             int nlat, int nlon, int leap,
             int *time, double *mld,
             double *taux, double *tauy ){
  int nmonth, retval, ncID, dataID, timeID, timeDimID, nn, mm, natts[2] = {0,0};
  size_t pos[2]={0,0}, start[3] = {0, (size_t) nlat, (size_t) nlon}, count[3]={1,1,1};
  size_t nt;
  int ntcum = 0;
  double *mmld, offset = 0.0, factor = 1.0;
  char filepath[MAXCHARLEN], timeDimName[NC_MAX_NAME + 1], attname[MAXCHARLEN];
  time_t itime[14];
  struct tm tcon;

  if (DBGFLG>2) { printf("  getdata: set MLD time\n"); fflush(NULL);}

  // set basic time constructor
  tcon.tm_year = (int) params[10] - 1900;
  tcon.tm_mon = 0;
  tcon.tm_mday = 0;
  tcon.tm_wday = 0;
  tcon.tm_yday = 0;
  tcon.tm_hour = 0;
  tcon.tm_min = 0;
  tcon.tm_sec = 0;
  // tcon.tm_gmtoff = 0;
  tcon.tm_isdst = 0;

  // take into account the reverse ordering in the mixed layer data
  pos[0]  = (size_t) nlat;
  pos[1]  = (size_t) nlon;

  mmld = dalloc(mmld, 16);

  // set hours of mid months
  tcon.tm_year = (int) params[10] - 1901;
  tcon.tm_mon = 10;
  tcon.tm_mday = 15;
  itime[0] = mktime(&tcon);

  tcon.tm_mon = 11;
  itime[1] = mktime(&tcon);

  tcon.tm_year = (int) params[10] - 1900;
  for (nn=2; nn<14; nn++){
    tcon.tm_mon = (nn - 1) % 12;
    tcon.tm_mday = 15;

    itime[nn] = mktime(&tcon);
  }

  tcon.tm_year = (int) params[10] - 1899;
  tcon.tm_mon = 0;
  tcon.tm_mday = 15;
  itime[14] = mktime(&tcon);

  tcon.tm_mon = 1;
  itime[15] = mktime(&tcon);

  if (DBGFLG>2) { printf("  getdata: load MLD\n"); fflush(NULL);}
  // load mixed layer depth
  for (nmonth=0; nmonth<12; nmonth++){
    sprintf(filepath, paths[2], nmonth+1);
    // open regridded netCDF file
    if ((retval = nc_open(filepath, NC_NOWRITE, &ncID))) ERR(retval);
    // get id of dataset
    if ((retval = nc_inq_varid(ncID, "DEPTH_MIXED_LAYER", &dataID))) ERR(retval);
    // get variable at position (nlon, nlat) and write into array
    if ((retval = nc_get_var1_double(ncID, dataID, pos, &mmld[nmonth+2]))) ERR(retval);
    // close nc file
    if ((retval = nc_close(ncID))) ERR(retval);
  }
  mmld[0]   = mmld[12]; // set november of previous year
  mmld[1]   = mmld[13]; // set december of previous year
  mmld[14]  = mmld[2]; // set january of next year
  mmld[15]  = mmld[3]; // set february of next year

  if (DBGFLG>2) { printf("  getdata: interpolate MLD in time\n"); fflush(NULL);}
  for (nn=0;nn<10248+leap*24;nn++){
    // check position relative to mid month times
    for (mm=1;mm<15;mm++) {
      if (itime[mm]>time[nn]){
        break;
      }
    }
    // interpolate
    mld[nn]   = (mmld[mm-1] +
                 (mmld[mm]-mmld[mm-1])/
                 (itime[mm]-itime[mm-1])*
                 (time[nn]-itime[mm-1]));
  }
  
  // load wind stress data
  if (DBGFLG>2) { printf("  getdata: load stress data\n"); fflush(NULL);}

    // iterate over months

  nt = 0;
  for (nmonth=0; nmonth<14; nmonth++) {
      if (nmonth==0) {
          sprintf(filepath, paths[1], (int) params[10]-1, 12);
      } else if (nmonth==13){
          sprintf(filepath, paths[1], (int) params[10]+1, 1);
      } else {
          sprintf(filepath, paths[1], (int) params[10], nmonth);
      }

      // open netcdf file
    if ((retval = nc_open(filepath, NC_NOWRITE, &ncID))) ERR(retval);
    // get time dimension
    if ((retval = nc_inq_varid(ncID, TIME, &timeID))) ERR(retval);
    if ((retval = nc_inq_vardimid (ncID, timeID, &timeDimID))) ERR(retval);
    if ((retval = nc_inq_dim(ncID, timeDimID, timeDimName, &nt))) ERR(retval);

    count[0] = nt;
    double *aux;
    aux = dalloc(aux, nt);
    for (nn=0;nn<nt;nn++){
      aux[nn]=0.0;
    }

    // get id of x dataset
    if ((retval = nc_inq_varid(ncID, "uflx", &dataID))) ERR(retval);
    // get variable at position (nlon, nlat) and write into array
    if ((retval = nc_get_vara_double(ncID, dataID, start, count, &aux[0]))) ERR(retval);


    // check how many attributes it has
    if ((retval = nc_inq_varnatts	(ncID, dataID, &natts[0]))) ERR(retval);
    // iterate over all attributes
    for (nn=0; nn<natts[0]; nn++){
      // read attribute name
      if ((retval = nc_inq_attname(ncID, dataID, nn, &attname[0]))) ERR(retval);
      // read offset if exists
      if (strcmp(attname, "add_offset")==0){
        if ((retval = nc_get_att_double(ncID, dataID, attname, &offset))) ERR(retval);
      }
      // read factor if exists
      if (strcmp(attname, "scale_factor")==0){
        if ((retval = nc_get_att_double(ncID, dataID, attname, &factor))) ERR(retval);
      }
    }

    if ((retval = nc_close(ncID))) ERR(retval);

    // set values according to scaling and offsets
    for (nn = 0; nn<nt; nn++) {
      taux[ntcum + nn] = factor * aux[nn] + offset;
    }

    // reset factor and offset
    offset = 0.0;
    factor = 1.0;

    if ((retval = nc_open(filepath, NC_NOWRITE, &ncID))) ERR(retval);

    // get id of y dataset
    if ((retval = nc_inq_varid(ncID, "vflx", &dataID))) ERR(retval);
    // get variable at position (nlon, nlat) and write into array
    if ((retval = nc_get_vara_double(ncID, dataID, start, count, &aux[0]))) ERR(retval);
    // check how many attributes it has
    if ((retval = nc_inq_varnatts	(ncID, dataID, &natts[0]))) ERR(retval);
    // iterate over all attributes
    for (nn=0; nn<natts[0]; nn++) {
      // read attribute name
      if ((retval = nc_inq_attname(ncID, dataID, nn, &attname[0]))) ERR(retval);
      // read offset if exists
      if (strcmp(attname, "add_offset") == 0) {
        if ((retval = nc_get_att_double(ncID, dataID, attname, &offset))) ERR(retval);
      }
      // read factor if exists
      if (strcmp(attname, "scale_factor") == 0) {
        if ((retval = nc_get_att_double(ncID, dataID, attname, &factor))) ERR(retval);
      }
    }

    // close nc file
    if ((retval = nc_close(ncID))) ERR(retval);

    // set values according to scaling and offsets
    for (nn = 0; nn < nt; nn++) {
      tauy[ntcum + nn] = factor * aux[nn] + offset;
    }

    ntcum += (int) nt;

    // free auxiliaries
    free(aux);
  }

  // correct stress values
  if ((int) params[4]!=0) {
    if (DBGFLG>2) { printf("  getdata: correcting NCEP-CFSR stress data\n"); fflush(NULL);}
    correct(taux, ntcum);
    correct(tauy, ntcum);
  }

  free(mmld);

  if (DBGFLG>2) { printf("  getdata: return to main\n"); fflush(NULL);}
}

dat2d_2 initdamping(dat2d *lsmask){
  int nn, mm;
  int ncID, varID[3], dimID, retval;
  size_t dampNLat;
  dat2d_2 data;
  char filepath[] = DAMPPATH, dimName[NC_MAX_NAME];
  double *dampLat, *dampWorld, *dampNA, *rr, *rrNA;

  if (DBGFLG>2) { printf("  initdamping: reading damping time scales\n"); fflush(NULL);}

  // open netcdf file
  if ((retval = nc_open(filepath, NC_NOWRITE, &ncID))) ERR(retval);

  // find variable IDs
  if ((retval = nc_inq_varid(ncID, LATS, &varID[0]))) ERR(retval);
  if ((retval = nc_inq_varid(ncID, DAMP, &varID[1]))) ERR(retval);
  if ((retval = nc_inq_varid(ncID, DAMPNA, &varID[2]))) ERR(retval);

  // find dimension bounds
  if ((retval = nc_inq_vardimid(ncID, varID[0], &dimID))) ERR(retval);
  if ((retval = nc_inq_dim(ncID, dimID, dimName, &dampNLat))) ERR(retval);

  // allocate meomory for temporary arrays
  dampLat = dalloc(dampLat, dampNLat);
  dampWorld = dalloc(dampWorld, dampNLat);
  dampNA = dalloc(dampNA, dampNLat);

  rr = dalloc(rr, (size_t) lsmask->nlat);
  rrNA = dalloc(rrNA, (size_t) lsmask->nlat);

  // read data
  if ((retval = nc_get_var_double(ncID, varID[0], &dampLat[0]))) ERR(retval);
  if ((retval = nc_get_var_double(ncID, varID[1], &dampWorld[0]))) ERR(retval);
  if ((retval = nc_get_var_double(ncID, varID[2], &dampNA[0]))) ERR(retval);

  // close nc file
  if ((retval = nc_close(ncID))) ERR(retval);

  // interpolate to model grid
  for (nn=0; nn<lsmask->nlat; nn++) {
    // check position relative to given data latitude
    for (mm = 1; mm < dampNLat; mm++) {
      if (dampLat[mm] > lsmask->lat[nn]) {
        break;
      }
    }
    // inter- / extrapolate
    if (lsmask->lat[nn] < dampLat[0]){ // extrapolate constant in south
      rr[nn] = dampWorld[0];
      rrNA[nn] = dampNA[0];
    } else if (lsmask->lat[nn] > dampLat[dampNLat-1]) { // extrapolate constant in north
      rr[nn] = dampWorld[dampNLat-1];
      rrNA[nn] = dampNA[dampNLat-1];
    } else { // if north, interpolate
      rr[nn] = (dampWorld[mm - 1] +
              (dampWorld[mm] - dampWorld[mm - 1]) /
              (dampLat[mm] - dampLat[mm - 1]) *
              (lsmask->lat[nn] - dampLat[mm - 1]));
      rrNA[nn] = (dampNA[mm - 1] +
              (dampNA[mm] - dampNA[mm - 1]) /
              (dampLat[mm] - dampLat[mm - 1]) *
              (lsmask->lat[nn] - dampLat[mm - 1]));

    }
  }

  // free temporary arrays
  free(dampLat);
  free(dampWorld);
  free(dampNA);

  // assign data to struct
  data.nx = lsmask->nlat;
  data.xx = lsmask->lat;
  data.y1 = rr;
  data.y2 = rrNA;

  return(data);
}

void getdataHybrid(double* params, char **paths,
                   dat2d *lsmask, dat1d *lh, dat1d *ww, dat1d *NN,
                   int ny, int nymin, int nx, int nxmin, int leap, int nlat5, int slat5){
  // indices
  size_t days[12] = {31, 28 + (size_t) leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int nmonth, nt, nint;
  // nc variables
  int retval, ncID, varID[2];
  size_t start[3], count[3], pos[2];
  char filepath[MAXCHARLEN];
  // time variables
  struct tm tcon;

  // initialize auxiliary variables
  double *aux;
  int *buoyancyTime, startIndex=0;
  pos[0] = (size_t) ny;
  pos[1] = (size_t) nx;
  start[1] = (size_t) ny;
  start[2] = (size_t) nx;
  start[0] = 0;
  count[1] = count[2] = 1;

  // set lh time for interpolation
  buoyancyTime = ialloc(buoyancyTime, 14);
  for (nt=0; nt<14; nt++){
    if (nt==0){
      tcon.tm_year = (int) params[10] - 1901;
      tcon.tm_mon = 11;
    } else if (nt==13) {
      tcon.tm_year = (int) params[10] - 1899;
      tcon.tm_mon = 0;
    } else {
      tcon.tm_year = (int) params[10] - 1900;
      tcon.tm_mon = nt-1;
    }
    tcon.tm_mday = 15;
    tcon.tm_hour = 1;
    tcon.tm_min = 0;
    tcon.tm_sec = 0;
    // tcon.tm_gmtoff = 0;
    tcon.tm_isdst = 0;
    buoyancyTime[nt] =  (int) mktime(&tcon);
  }

  if (DBGFLG>2) { printf("  getdataHybrid: load buoyancy frequency data\n"); fflush(NULL);}
  // load N data
  aux = dalloc(aux, 14);
  for (nmonth=1; nmonth<13; nmonth++) {
    // init aux at nmonth
    aux[nmonth]=0.0;

    // set filepath
    sprintf(filepath, paths[0], nmonth);

    // open file
    if ((retval = nc_open(filepath, NC_NOWRITE, &ncID))) ERR(retval);

    // get varID
    if ((retval = nc_inq_varid(ncID, BUOYANCY, &varID[0]))) ERR(retval);

    // get variable
    if ((retval = nc_get_var1_double(ncID, varID[0], pos, &aux[nmonth]))) ERR(retval);

    // close file
    if ((retval = nc_close(ncID))) ERR(retval);
  }
  // set previous December and following January
  aux[0] = aux[12];
  aux[13] = aux[1];

  if (DBGFLG>2) { printf("  getdataHybrid: interpolate buoyancy frequency data\n"); fflush(NULL);}
  // interpolate buoyancy frequency in time
  for (nt=0; nt<lh->ntime; nt++) {
    // check position relative to monthly values
    for (nint = 1; nint < 14; nint++) {
      if (buoyancyTime[nint] > lh->time[nt]) {
        break;
      }
    }
    // linearly interpolate and set to data array
    NN->data[nt] = (aux[nint - 1] +
                   (aux[nint] - aux[nint - 1]) /
                   (buoyancyTime[nint] - buoyancyTime[nint - 1]) *
                   (lh->time[nt] - buoyancyTime[nint - 1]));
  }

  if (DBGFLG>2) { printf("  getdataHybrid: load vertical velocities and horizontal length scales\n"); fflush(NULL);}
  // load vertical velocity and horizontal length scale
  for (nmonth=0; nmonth<12; nmonth++){
    startIndex = sum(days, nmonth) * 24;
    count[0] = days[nmonth]*24;

    // set filepath to datafiles written earlier

    if ((nlat5==NC_FILL_INT)|(slat5==NC_FILL_INT)) {
      if (lsmask->lat[ny] > 0) {
        sprintf(filepath, paths[5], (int) params[10], nmonth+1);
      } else if (lsmask->lat[ny] < 0) {
        sprintf(filepath, paths[6], (int) params[10], nmonth+1);
      }
      start[1] = (size_t) (ny-nymin);
    } else {
      if (lsmask->lat[ny] > lsmask->lat[nlat5]) {
        sprintf(filepath, paths[5], (int) params[10], nmonth + 1);
        if (lsmask->lat[nymin] > lsmask->lat[nlat5]) {
          start[1] = (size_t) (ny - nymin);
        } else if (lsmask->lat[nymin] <= lsmask->lat[nlat5]) {
          start[1] = (size_t) (ny - nlat5);
        }
      } else if (lsmask->lat[ny] < lsmask->lat[slat5]) {
        sprintf(filepath, paths[6], (int) params[10], nmonth + 1);
        if (lsmask->lat[nymin] < lsmask->lat[slat5]) {
          start[1] = (size_t) (ny - nymin);
        } else if (lsmask->lat[nymin] >= lsmask->lat[slat5]) {
          start[1] = (size_t) (ny - slat5);
        }
      } else {
        // if lat of point is in forbidden latitude band throw an error
        GENERR
      }
    }
    start[2] = (size_t) (nx - nxmin);

    // open file
    if ((retval = nc_open(filepath, NC_NOWRITE, &ncID))) ERR(retval);

    // get varID
    if ((retval = nc_inq_varid(ncID, ZVEL, &varID[0]))) ERR(retval);
    if ((retval = nc_inq_varid(ncID, LH, &varID[1]))) ERR(retval);

    // get variable
    if ((retval = nc_get_vara_double(ncID, varID[0], start, count, &ww->data[startIndex]))) ERR(retval);
    if ((int) params[3]==1) {
      if ((retval = nc_get_vara_double(ncID, varID[1], start, count, &lh->data[startIndex]))) ERR(retval);
    }

    // close file
    if ((retval = nc_close(ncID))) ERR(retval);
  }

  if ((int) params[3]!=1) {
    for (nt=0; nt<lh->ntime; nt++){
      lh->data[nt] = params[16];
    }
  }

  free(aux);
  free(buoyancyTime);
}