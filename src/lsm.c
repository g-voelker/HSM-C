#include <stdio.h>
#include "../lib/complex.h"
#include "../lib/alloc_space.h"
#include "header.h"
#include "input.h"
#include <netcdf.h>
#include <stdlib.h>

struct dat {
  int nlat, nlon;
  double *lat, *lon;
  double **data;
};

struct dat lsm(char *filename) {

  // initialize variables
  int dataID, latID, lonID, ncID;
  size_t NLAT, NLON;
  double *lat, *lon, ***data;

  // error handling, dimensions
  int retval, latDimID[1], lonDimID[1];
  char latDimName[NC_MAX_NAME + 1], lonDimName[NC_MAX_NAME + 1];
  struct dat lsm;

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

  int datadimID[3], nn, ndim[1];
  size_t datadim;
  char datadimName[NC_MAX_NAME + 1];
  if ((retval = nc_inq_vardimid (ncID, dataID, datadimID))) ERR(retval);

//  for (nn=0; nn<=2; nn++) {
//    if ((retval = nc_inq_dim(ncID, datadimID[nn], datadimName, &datadim))) ERR(retval);
//    printf(datadimName);printf(": ");fflush(NULL);
//    printf("%d\n", datadim);
//    fflush(NULL);
//  }

  if (DBGFLG>2) { printf("lsm: allocate space\n"); fflush(NULL);}
  
  // allocate space for vectors
  lat = dvector(0, NLAT - 1);
  lon = dvector(0, NLON - 1);
  data = dmatrix3(0, 0, 0, NLAT - 1, 0, NLON - 1);
  free_dvector(lat, 0, NLAT - 1);
  free_dvector(lon, 0, NLON - 1);
  free_dmatrix3(data, 0, 0, 0, NLAT - 1, 0, NLON - 1);

  printf("Das wars schon");

  lat = dvector(0, NLAT - 1);
  lon = dvector(0, NLON - 1);
  data = dmatrix3(0, 0, 0, NLAT - 1, 0, NLON - 1);

//  double data[1][NLAT][NLON];

  if (DBGFLG>2) { printf("lsm: read land sea mask\n"); fflush(NULL);}

  size_t vstart[] = {0, 0, 0};
  size_t vlen[] = {1, NLAT - 1, NLON - 1};
  ptrdiff_t stride[] = {1, 1, 1};

  // read data arrays
  if ((retval = nc_get_var_double(ncID, latID, &lat[0]))) ERR(retval);
  if ((retval = nc_get_var_double(ncID, lonID, &lon[0]))) ERR(retval);
  if ((retval = nc_get_vars_double(ncID, dataID, vstart, vlen, stride, &data[0][0][0]))) ERR(retval);

  if (DBGFLG>2) { printf("lsm: close netcdf file\n"); fflush(NULL);}

  // close netcdf file
  if ((retval = nc_close(ncID))) ERR(retval);

  lsm.nlat = NLAT;
  lsm.nlon = NLON;
  lsm.lat = lat;
  lsm.lon = lon;
  lsm.data = data[0];

//  free_dvector(lat, 0, NLAT - 1);
//  free_dvector(lon, 0, NLON - 1);
//  free_dmatrix3(data, 0, 0, 0, NLAT - 1, 0, NLON - 1);

  return lsm;
}

