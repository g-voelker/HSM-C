
struct dat1D {
  int ntime;
  double *time;
  double *data;
};

struct dat2D {
  int nlat, nlon;
  double *lat, *lon;
  double **data;
};

struct dat3D {
  int nt, nlat, nlon;
  double *lat, *lon, *time;
  double ***data;
};

typedef struct dat1D dat1d;
typedef struct dat2D dat2d;
typedef struct dat3D dat3d;