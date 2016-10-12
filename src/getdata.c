#include <grib_api.h>
#include <netcdf.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "input.h"
#include "header.h"

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
  
  double temp[] = input;
  int nn;
  for (nn=1; nn<length; nn++) {
    input[nn] = temp[n]*(nn%6+1) - temp[nn-1]*(nn%6);
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
  int nmonth, retval, pos[2], ncID, dataID, leap, nn;
  double *mmld, *mtime;
  
  pos[1]  = nlat;
  pos[0]  = nlon;
  mmld     = dvector(0,13)
  mtime   = dvector(0,13)
  
  // load mixed layer depth
  for (nmonth=1; nmonth<=12; nmonth++){
    char filepath[] = MLDPATH % nmonth;
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
  
  // interpolate mixed layer depth in time
  // set hours of mid months
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
  
  for (nn=0;nn<=8760+leap*24;n++){
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
