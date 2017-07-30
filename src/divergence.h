//
// Created by georg on 24/11/16.
//

extern double dist(double lon1, double lon2, double lat1, double lat2);
extern void divergence(double *params, char **paths, dat2d *lsmask,
                       int nxmin, int nxmax, int nymin, int nymax, int leap, int hemflag);