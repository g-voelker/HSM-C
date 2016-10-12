

#include <math.h>
#include <stdio.h>
#include "complex.h"
#include "alloc_space.h"
#include "header.h"
#include <netcdf.h>
#include <stdlib.h>

void lsm(char *filename, int ncID, int dataID, int latID, int lonID,
		 size_t NLAT, size_t NLON,
         double ***data, double *lat, double *lon) {
  // error handling, dimensions
  int retval, latDimID[1], lonDimID[1], nn, mm, nt;
  char latDimName[NC_MAX_NAME], lonDimName[NC_MAX_NAME];
  
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
  lat = dvector(0, NLAT-1);
  lon = dvector(0, NLON-1);
  data = dmatrix3(0, 0, 0, NLAT-1, 0, NLON-1);

  if (DBGFLG>2) { printf("lsm: read land sea mask\n"); fflush(NULL);}
  
  // read data arrays
  if ((retval = nc_get_var_double(ncID, latID, &lat[0]))) ERR(retval);
  if ((retval = nc_get_var_double(ncID, lonID, &lon[0]))) ERR(retval);
  if ((retval = nc_get_var_double(ncID, dataID, &data[0][0][0]))) ERR(retval);
  
  // close netcdf file
  if ((retval = nc_close(ncID))) ERR(retval);
}

