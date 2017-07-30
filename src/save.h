//
// Created by georg on 21/11/16.
//

extern int sum(size_t *array, int num);

extern void initnc(double *params, char **paths, dat2d *lsmask,
                   int* time, int xmin, int xmax, int ymin, int ymax, int leap, int *nlat5, int *slat5);

extern void savePoint(double *params, char **paths, dat2d *lsmask,
                      double* uu, double* vv, double* mld,  double* taux, double* tauy,
                      int* time, int nxmin, int nx, int nymin, int nymax, int ny, int leap, int nlat5, int slat5);

extern void savelh(double *params, char **paths, dat2d *lsmask,
                   double ***lh, int *time,
                   int nxmin, int nxmax, int nymin, int nymax, int leap, int hemflag);

extern void savePointHybrid(double *params, char **paths,
                            dat2d *lsmask, dat1d *Eout,
                            int ny, int nymin, int nymax, int nx, int nxmin, int leap, int nlat5, int slat5);