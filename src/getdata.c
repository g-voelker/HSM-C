#include <grib_api.h>
#include <netcdf.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "input.h"
#include "header.h"
#include "../lib/complex.h"
#include "../lib/alloc_space.h"
#include <grib_api.h>
//#include <unistd.h> (sleep)


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
  
  double *temp = input;
  int nn;
  for (nn=1; nn<length; nn++) {
    input[nn] = temp[nn]*(nn%6+1) - temp[nn-1]*(nn%6);
  }
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
             double *time, double *mld,
             double *taux, double *tauy ){
  int nmonth, retval, ncID, dataID, timeID, timeDimID, nn, mm, natts[2] = {0,0};
  size_t pos[2]={0,0}, start[3] = {0, (size_t) nlat, (size_t) nlon}, count[3]={1,1,1};
  size_t nt;
  int ntcum = 0;
  double *mmld, *mtime, offset = 0.0, factor = 1.0;
  char filepath[MAXCHARLEN], timeDimName[NC_MAX_NAME + 1], attname[MAXCHARLEN];
  time_t itime[14];
  struct tm tcon;
//  FILE* in = NULL;
//  grib_handle *handle = NULL;
//  char shortName[MAXCHARLEN];
//  int grib_count=0;
//  size_t vlen=MAXCHARLEN;
//  long Ni=0, Nj=0;

  if (DBGFLG>2) { printf("getdata: set MLD time\n"); fflush(NULL);}

  // set basic time constructor
  tcon.tm_year = YEAR - 1900;
  tcon.tm_mon = 0;
  tcon.tm_mday = 0;
  tcon.tm_wday = 0;
  tcon.tm_yday = 0;
  tcon.tm_hour = 0;
  tcon.tm_min = 0;
  tcon.tm_sec = 0;
  tcon.tm_gmtoff = 0;
  tcon.tm_isdst = 0;

  pos[0]  = (size_t) nlat;
  pos[1]  = (size_t) nlon;

  mmld     = dvector(0,13);

  // set hours of mid months
  tcon.tm_year = YEAR - 1901;
  tcon.tm_mon = 11;
  tcon.tm_mday = 15;
  itime[0] = mktime(&tcon);

  tcon.tm_year = YEAR - 1900;
  for (nn=1; nn<13; nn++){
    tcon.tm_mon = (nn - 1) % 12;
    tcon.tm_mday = 15;

    itime[nn] = mktime(&tcon);
  }

  tcon.tm_year = YEAR - 1899;
  tcon.tm_mon = 0;
  tcon.tm_mday = 15;
  itime[13] = mktime(&tcon);

//  for (nn=0; nn<14; nn++){
//    printf("%s\n", ctime(&itime[nn]));fflush(NULL);
//    printf("%d\n", (int) (itime[nn] - itime[0]));fflush(NULL);
//  }

  if (DBGFLG>2) { printf("getdata: load MLD\n"); fflush(NULL);}
  // load mixed layer depth
  for (nmonth=1; nmonth<=12; nmonth++){
    sprintf(filepath, MLDPATH, nmonth);
    // open regridded netCDF file
    if ((retval = nc_open(filepath, NC_NOWRITE, &ncID))) ERR(retval);
    // get id of dataset
    if ((retval = nc_inq_varid(ncID, "DEPTH_MIXED_LAYER", &dataID))) ERR(retval);
    // get variable at position (nlon, nlat) and write into array
    if ((retval = nc_get_var1_double(ncID, dataID, pos, &mmld[nmonth]))) ERR(retval);
    // close nc file
    if ((retval = nc_close(ncID))) ERR(retval);
  }
  mmld[0]   = mmld[12];
  mmld[13]  = mmld[1];

  for (nn=0;nn<=8760+leap*24;nn++){
    // check position relative to mid month times
    for (mm=1;mm<=12;mm++) {
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
  if (DBGFLG>2) { printf("getdata: load stress data\n"); fflush(NULL);}

  // iterate over months

  nt = 0;
  for (nmonth=1; nmonth<=12; nmonth++) { // this MUST start at 1.
    sprintf(filepath, STRSPATH, YEAR, nmonth);
    // open netcdf file
    if ((retval = nc_open(filepath, NC_NOWRITE, &ncID))) ERR(retval);
    // get time dimension
    if ((retval = nc_inq_varid(ncID, "time", &timeID))) ERR(retval);
    if ((retval = nc_inq_vardimid (ncID, timeID, &timeDimID))) ERR(retval);
    if ((retval = nc_inq_dim(ncID, timeDimID, timeDimName, &nt))) ERR(retval);

    count[0] = nt;
    double aux[nt];
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
  }

  // correct stress values
  if (STRSCOR!=0) {
    if (DBGFLG>2) { printf("getdata: correcting NCEP-CFSR stress data\n"); fflush(NULL);}
    correct(taux, ntcum);
    correct(tauy, ntcum);
  }

  if (DBGFLG>2) { printf("getdata: return to main\n"); fflush(NULL);}
}
