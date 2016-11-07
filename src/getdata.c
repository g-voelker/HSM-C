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
  free(temp);
}

/* data to get, given nlat, nlon:
 * mld (in m)
 * tau_x (in N/m^2)
 * tau_y (in N/m^2)
 * time axis (in seconds from 1/1/YEAR 00:00:00)
 * 
 * modifications that need to be done:
 * interpolate mld to given year's timeline
 * correct stress values as they come in averages
 */

void getdata(int nlat, int nlon, int leap,
             double *time, double *mld,
             double *taux, double *tauy ){
  int nmonth, retval, ncID, dataID, nn, mm;
  size_t pos[2];
  double *mmld, *mtime;
  char filepath[MAXCHARLEN];
  time_t itime[14];
  struct tm tcon;

  if (DBGFLG>2) { printf("getdata: set time constructor\n"); fflush(NULL);}
  // set basic time constructor
  tcon.tm_year = YEAR - 1900;
  tcon.tm_mon = 0;
  tcon.tm_mday = 0;
  tcon.tm_hour = 0;
  tcon.tm_min = 0;
  tcon.tm_sec = 0;
  
  pos[1]  = (size_t) nlat;
  pos[0]  = (size_t) nlon;
  mmld     = dvector(0,13);
  mtime   = dvector(0,13);

  if (DBGFLG>2) { printf("getdata: set MLD time\n"); fflush(NULL);}
  // set hours of mid months
  tcon.tm_year = YEAR - 1901;
  tcon.tm_mon = 11;
  tcon.tm_mday = 15;
  itime[0] = mktime(&tcon);

  if (DBGFLG>2) { printf("getdata: break point 1\n"); fflush(NULL);}

  tcon.tm_year = YEAR - 1900;
  for (nn=1; nn<13; nn++){
    tcon.tm_mon = (nn - 1) % 12;
    tcon.tm_mday = 15;

    itime[nn] = mktime(&tcon);
  }

  if (DBGFLG>2) { printf("getdata: break point 2\n"); fflush(NULL);}

  tcon.tm_year = YEAR - 1899;
  tcon.tm_mon = 0;
  tcon.tm_mday = 15;
  itime[13] = mktime(&tcon);

  if (DBGFLG>2) { printf("getdata: break point 3\n"); fflush(NULL);}

  for (nn=0; nn<14; nn++){
    printf("%s\n", ctime(&itime[nn]));fflush(NULL);
  }

  if (DBGFLG>2) { printf("getdata: break point 4\n"); fflush(NULL);}

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
  

  mtime[0]  = -31/2*24;
  mtime[1]  = 31/2*24;
  mtime[2]  = mtime[1]*2 + (28 + leap)/2*24;
  mtime[3]  = (59 + leap)*24 + 31/2*24;
  mtime[4]  = (90 + leap)*24 + 30/2*24;
  mtime[5]  = (120 + leap)*24 + 31/2*24;
  mtime[6]  = (151 + leap)*24 + 30/2*24;
  mtime[7]  = (181 + leap)*24 + 31/2*24;
  mtime[8]  = (212 + leap)*24 + 31/2*24;
  mtime[9]  = (243 + leap)*24 + 30/2*24;
  mtime[10] = (273 + leap)*24 + 31/2*24;
  mtime[11] = (304 + leap)*24 + 30/2*24;
  mtime[12] = (334 + leap)*24 + 31/2*24;
  mtime[13] = (365 + leap)*24 + 31/2*24;
  
  for (nn=0;nn<=8760+leap*24;nn++){
    // set time entry (iterates over hours)
    time[nn]  = nn;
    
    // check position relative to mid month times
    for (mm=1;mm<=12;mm++) {
      if (mtime[mm]>time[nn]){
        break;
      }
    }
    
    // interpolate
    mld[nn]   = (mmld[mm-1] +
                 (mmld[mm]-mmld[mm-1])/
                 (mtime[mm]-mtime[mm-1])*
                 (time[nn]-mtime[mm-1]));
  }
  
  // load wind stress data
  
  
  
  
}
