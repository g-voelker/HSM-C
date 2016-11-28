#include "../lib/constants.h"

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define DIMERR(e) {printf("Error: too few points in given domain. Hybrid model aborted.\n"); exit(e);}
#define DEG2RAD(deg) ((deg)*PI/180.0)
#define RAD2DEG(rad) {(rad)*180.0/PI);}
#define MAXCHARLEN 100

// the length of the DFTs is set to 2^5 * 7^3 in accordance with the maximum cutoff for long damping times.
// this value is set to preserve performance while making sure the cutoff is long enough
// edit only if you know what you do.
#define DFFT_LEN 10976
#define CHUNK_LAT 2
#define CHUNK_LON 2

// macros for variable / attribute name space
#define UNITS "units"
#define LONGNAME "long_name"

#define TIME "time"
#define LATS "latitude"
#define LONS "longitude"
#define XVEL "u"
#define XVEL_LONG "x-component of velocity in mixed layer"
#define YVEL "v"
#define YVEL_LONG "y-component of velocity in mixed layer"
#define ZVEL "w"
#define ZVEL_LONG "vertical velocity component at base of  mixed layer"
#define MLD "mld"
#define MLD_LONG "mixed layer depth"

#define DEGREES_NORTH "degrees_north"
#define DEGREES_EAST "degrees_east"
#define HOURS "hours since 1900-01-01 00:00:0.0"
#define MPS "meter/second"
#define METER "meter"