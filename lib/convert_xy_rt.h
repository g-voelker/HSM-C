extern void convert_xy_rt(
   double **fxy, int ox, int nx, int oy, int ny,
   double xmin, double xmax, double ymin, double ymax,
   double **frt, int or, int nr, int ot, int nt,
   double rmin, double rmax, double tmin, double tmax);
extern void convert_rt_xy(
   double **frt, int or, int nr, int ot, int nt,
   double rmin, double rmax, double tmin, double tmax,
   double **fxy, int ox, int nx, int oy, int ny,
   double xmin, double xmax, double ymin, double ymax);
