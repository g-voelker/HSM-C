#include "../lib/constants.h"

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
#define DIMERR(e) {printf("Error: too few points in given domain. Hybrid model aborted.\n"); exit(e);}
#define DEG2RAD(deg) ((deg)*PI/180.0)
#define RAD2DEG(rad) {(rad)*180.0/PI);}
#define MAXCHARLEN 100

// the length of the DFTs is set to 2^5 * 7^3 in accordance with the maximum cutoff for long damping times.
// this value is set to preserve performance while making sure the cutoff is long enough
// edit only if you know what you do; values are related to ringing effects
#define DFFT_LEN 10976

// chunking settings for nc files saved
// there is a performance / diec space trade off
// at the moment it is optimized to performance
// edit only if you know what you do; too high values may lead to crashes
#define CHUNK_LAT 1
#define CHUNK_LON 1

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
#define TAUX "taux"
#define TAUX_LONG "Momentum flux, u component"
#define TAUY "tauy"
#define TAUY_LONG "Momentum flux, v component"
#define EIN "wind_work"
#define EIN_LONG "work done by wind on ocean surface"
#define EOUT "phi_IW"
#define EOUT_LONG "energy flux due to radiated IWs"

#define DEGREES_NORTH "degrees_north"
#define DEGREES_EAST "degrees_east"
#define HOURS "hours since 1900-01-01 00:00:0.0"
#define MPS "m/s"
#define METER "m"
#define NPM2 "N / m**2"
#define WPM2 "W / m**2"