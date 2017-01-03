//
// Created by georg on 21/11/16.
//

extern int sum(size_t *array, int num);

extern void initnc(dat2d *lsmask, int* time, int xmin, int xmax, int ymin, int ymax, int leap);

extern void savePoint(double* uu, double* vv, double* mld,  double* taux, double* tauy,
                      int* time, int nxmin, int nx, int nymin, int ny, int leap);

extern void savelh(double ***lh, int *time, int nxmin, int nxmax, int nymin, int nymax, int leap);

extern void savePointHybrid(dat1d *Eout, int ny, int nymin, int nx, int nxmin, int leap);