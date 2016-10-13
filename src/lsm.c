#include <stdio.h>
#include "../lib/complex.h"
#include "../lib/alloc_space.h"
#include "header.h"
#include "input.h"
#include <netcdf.h>
#include <stdlib.h>

struct dat2D {
  int nlat, nlon;
  double *lat, *lon;
  double **data;
};

struct dat2D lsm(char *filename) {

  // initialize variables
  int dataID, latID, lonID, ncID, nn, mm;
  size_t NLAT, NLON;
  double *lat, *latd, *lon, *data, **mask;

  // error handling, dimensions
  int retval, latDimID[1], lonDimID[1];
  char latDimName[NC_MAX_NAME + 1], lonDimName[NC_MAX_NAME + 1];
  struct dat2D lsm;

  if (DBGFLG>2) { printf("lsm: open netcdf file\n"); fflush(NULL);}

  // open file
  if ((retval = nc_open(filename, NC_NOWRITE, &ncID))) ERR(retval);
  
  if (DBGFLG>2) { printf("lsm: get variable ids\n"); fflush(NULL);}
  
  // determine variable IDs
  if ((retval = nc_inq_varid(ncID, "lat", &latID))) ERR(retval);
  if ((retval = nc_inq_varid(ncID, "lon", &lonID))) ERR(retval);
  if ((retval = nc_inq_varid(ncID, "LAND_L1", &dataID))) ERR(retval);
  
  if (DBGFLG>2) { printf("lsm: set variable bounds\n"); fflush(NULL);}
  
  // determine lengths of lat and lon
  if ((retval = nc_inq_vardimid (ncID, latID, latDimID))) ERR(retval);
  if ((retval = nc_inq_vardimid (ncID, lonID, lonDimID))) ERR(retval);
  if ((retval = nc_inq_dim (ncID, latDimID[0], latDimName, &NLAT))) ERR(retval);
  if ((retval = nc_inq_dim (ncID, lonDimID[0], lonDimName, &NLON))) ERR(retval);

  if (DBGFLG>2) { printf("lsm: allocate space\n"); fflush(NULL);}
  
  // allocate space for vectors
  latd = dvector(0, NLAT - 1);
  lat = dvector(0, NLAT - 1);
  lon = dvector(0, NLON - 1);
  data = dvector(0, NLAT*NLON - 1);
  mask = dmatrix2(0, NLON - 1, 0, NLAT - 1);

  if (DBGFLG>2) { printf("lsm: read land sea mask\n"); fflush(NULL);}

  // read data arrays
  if ((retval = nc_get_var_double(ncID, latID, &latd[0]))) ERR(retval);
  if ((retval = nc_get_var_double(ncID, lonID, &lon[0]))) ERR(retval);
  if ((retval = nc_get_var_double(ncID, dataID, &data[0]))) ERR(retval);



  for (nn=0; nn<NLAT; nn++){
    lat[NLAT - 1 - nn] = latd[nn];
    for (mm=0; mm<NLON; mm++) {
      mask[mm][NLAT - 1 - nn] = data[mm + nn * NLON];
    }
  }

  if (DBGFLG>2) { printf("lsm: close netcdf file\n"); fflush(NULL);}

  // close netcdf file
  if ((retval = nc_close(ncID))) ERR(retval);

  lsm.nlat = (int) NLAT;
  lsm.nlon = (int) NLON;
  lsm.lat = lat;
  lsm.lon = lon;
  lsm.data = mask;

  // pointers are assigned to struct so can't free them.
//  free_dvector(lat, 0, NLAT - 1);
//  free_dvector(lon, 0, NLON - 1);
//  free_dmatrix3(data, 0, 0, 0, NLAT - 1, 0, NLON - 1);

  return lsm;
}

