//
// Created by Georg Sebastian Voelker on 13/10/16.
//

extern void correct(double *input, int length);

extern void getdata(int nlat, int nlon, int leap,
                    int *time, double *mld,
                    double *taux, double *tauy );

extern dat2d_2 initdamping(dat2d *lsmask);

extern void getdataHybrid(dat2d *lsmask, dat1d *lh, dat1d *ww, dat1d *NN,
                          int ny, int nymin, int nx, int nxmin, int leap, int nlat5, int slat5);