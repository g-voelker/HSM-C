extern void writeXPlot(
   char *sfpw,
   double *x, double *f, int nc, int *npts);
extern void writeLogXPlot(
   char *sfpw,
   double *x, double *f, int nc, int *npts);
extern void writeLogLogXPlot(
   char *sfpw,
   double *x, double *f, int nc, int *npts);

extern void writeXYPlot(
   char *sfpw,
   double *x, double *y, double *data, int nc, int **nxypts, int *ctype);
extern void writeConPlot(
   char *sfpw,
   double **f, int nx, int ny, double *params);
extern void writeLogConPlot(
   char *sfpw,
   double **f, int nx, int ny, double *params);
