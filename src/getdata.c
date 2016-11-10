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
  int nmonth, retval, ncID, dataID, nn, mm, nt;
  size_t pos[2]={0,0};
  double *mmld, *mtime;
  char filepath[MAXCHARLEN];
  time_t itime[14];
  struct tm tcon;
  FILE* in = NULL;
  grib_handle *handle = NULL;
  char shortName[MAXCHARLEN];
  int grib_count=0;
  size_t vlen=MAXCHARLEN;
  long Ni=0, Nj=0;

  if (DBGFLG>2) { printf("getdata: set MLD time\n"); fflush(NULL);}

  // set basic time constructor
  tcon.tm_year = YEAR - 1900;
  tcon.tm_mon = 0;
  tcon.tm_mday = 0;
  tcon.tm_hour = 0;
  tcon.tm_min = 0;
  tcon.tm_sec = 0;

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
  // turn on multi message support
  grib_multi_support_on(0);

//  char value[MAXCHARLEN];
//  char* name_space="geography";
//  unsigned long key_iterator_filter_flags=GRIB_KEYS_ITERATOR_ALL_KEYS;

  // iterate over months
  nt = 0;
  for (nmonth=1; nmonth<=12; nmonth++) {
    sprintf(filepath, STRSPATH, YEAR, nmonth);
    in = fopen(filepath, "r");
    if (!in) {
      if (DBGFLG > 2) {
        printf("getdata: unable to open file %s\n", filepath);
        fflush(NULL);
      }
    }
    grib_count=0;
    while ((handle = grib_handle_new_from_file(0,in,&retval)) != NULL){
      if (grib_count==0) {
        GRIB_CHECK(grib_get_long(handle, "Ni", &Ni), 0);
        GRIB_CHECK(grib_get_long(handle, "Nj", &Nj), 0);
        printf("%d, %d\n", (int) Ni, (int) Nj);fflush(NULL);
      }

//      grib_keys_iterator* kiter=NULL;
//      kiter=grib_keys_iterator_new(handle,key_iterator_filter_flags,name_space);
//      printf("hour = %d\n", grib_count/2);

//      while(grib_keys_iterator_next(kiter)) {
//        const char *name = grib_keys_iterator_get_name(kiter);
//        if (strcmp(name, "iScansNegatively") == 0) {
//          vlen = MAXCHARLEN;
//          bzero(value, vlen);
//          GRIB_CHECK(grib_get_string(handle, name, value, &vlen), name);
//          printf("%s = %s\n", name, value);
//        }
//        printf("%s \n",name);
//      }

//      grib_keys_iterator_delete(kiter);


      GRIB_CHECK(grib_get_string(handle, "shortName", shortName, &vlen), "shortName");

      if (strcmp(shortName, "uflx") == 0){
        GRIB_CHECK(grib_get_double_element (handle, "values", (int) Ni*nlat + nlon, &taux[nt/2]), 0);
      } else if (strcmp(shortName, "vflx") == 0) {
        GRIB_CHECK(grib_get_double_element (handle, "values", (int) Ni*nlat + nlon, &tauy[nt/2]), 0);
      }

      grib_handle_delete(handle);
      grib_count++;
      nt++;
      if ((nt%20)==0) {
        printf("grib/total count: (%d, %d)\n", grib_count, nt);
      }
//      if (grib_count==2){break;}
    }
  }
  if (DBGFLG>2) { printf("getdata: return to main\n"); fflush(NULL);}
}
