extern double interpolate0(
   double xv, double x1, double x2, double y1, double y2);
extern double interpolate(
   double xv, 
   double x1, double x2, double x3, double y1, double y2,double y3);
extern void interpolate_list(
   double *x, double *y, int n, double xv, double *yv);

extern int findmax(
   double x1, double x2, double x3, double y1, double y2, double y3,
   double *xm, double *ym);

extern void myspline(
   double *x, double *y, int n, double yp1, double ypN, double *y2);
extern void mysplint(
   double *x, double *y, double *y2, int n, double xv, double *yv);


extern void interpolatefn(double *x,double *f,int n,
                          double *x2,double *f2,int n2);

