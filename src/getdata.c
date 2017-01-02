#include <netcdf.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../lib/dalloc.h"
#include "../lib/structs.h"
#include "input.h"
#include "header.h"
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

void getdata(int nlat, int nlon, int leap,
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
  tcon.tm_year = YEAR - 1900;
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
  tcon.tm_year = YEAR - 1901;
  tcon.tm_mon = 10;
  tcon.tm_mday = 15;
  itime[0] = mktime(&tcon);

  tcon.tm_mon = 11;
  itime[1] = mktime(&tcon);

  tcon.tm_year = YEAR - 1900;
  for (nn=2; nn<14; nn++){
    tcon.tm_mon = (nn - 1) % 12;
    tcon.tm_mday = 15;

    itime[nn] = mktime(&tcon);
  }

  tcon.tm_year = YEAR - 1899;
  tcon.tm_mon = 0;
  tcon.tm_mday = 15;
  itime[14] = mktime(&tcon);

  tcon.tm_mon = 0;
  itime[15] = mktime(&tcon);

  if (DBGFLG>2) { printf("  getdata: load MLD\n"); fflush(NULL);}
  // load mixed layer depth
  for (nmonth=0; nmonth<12; nmonth++){
    sprintf(filepath, MLDPATH, nmonth+1);
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
    if (DBGFLG>2) { printf("  nmonth: %d\n", nmonth); fflush(NULL);}

      if (nmonth==0) {
          sprintf(filepath, STRSPATH, YEAR-1, 12);
      } else if (nmonth==13){
          sprintf(filepath, STRSPATH, YEAR+1, 1);
      } else {
          sprintf(filepath, STRSPATH, YEAR, nmonth + 1);
      }
      if (DBGFLG>2) { printf(filepath); fflush(NULL);}

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
    if ((retval = nc_inq_varid(ncID, "p260062", &dataID))) ERR(retval);
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
    if ((retval = nc_inq_varid(ncID, "p260063", &dataID))) ERR(retval);
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
  if (STRSCOR!=0) {
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

void getlh(dat2d *lsmask, dat1d *lh, int nn, int mm){
  // get indicees of points to load

  // load matrix

  // analyse time step by time step

  // average time bins together

  // interpolate reduced time series

  // set values in struct

}