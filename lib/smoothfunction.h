extern void HeavysideAverage(
   double *xraw, double *fxraw, int nraw,
   double lx,double range,
   double *x, double *fx, int nn);

extern void GaussianAverage(
   double *xraw, double *fxraw, int nraw,
   double lx,double range,
   double *x, double *fx, double *sd, int nn);

extern void BestFitLineSmooth(
   double *xraw, double *fxraw, int nraw,
   double lx,double range,
   double *x, double *fx, int nn);

extern void BestFitLineFilter(
   double *xraw, double *fxraw, int nraw,
   double lx,double range, double nsdm, double nsdb,
   double *x, double *fx, int nn);

extern void BestFitLineDiff(
   double *xraw, double *fxraw, int nraw,
   double lx,double range,
   double *x, double *fx, int nn);

