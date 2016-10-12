extern void readXPlot(
   char *sfpr,
   double *x, double *fx, int *ncurves, int *npts, double *xyrng,
   int NC, int NN);
extern void readXPlot_noprompt(
   char *sfpr,
   double *x, double *fx, int *ncurves, int *npts, double *xyrng,
   int NC, int NN);
extern void readXYPlot(
   char *sfpr,
   double *x, double *y, double *data, int *ncurves, int **nxypts, 
   int NC, int NN);
extern void readXYPlot_noprompt(
   char *sfpr,
   double *x, double *y, double *data, int *ncurves, int **nxypts, 
   int NC, int NN);
extern void readConPlot(
   char *sfpr,
   double **fld, int *nxpts, int *nypts, double *params);
