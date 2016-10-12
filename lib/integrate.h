extern void *int_derivs(
   double xind, double y, double *yp, int npts,
   double *x, double *fx, double *fxspline);

extern void myrk4(
   double y, double dydx, double x, double h, double *yout,
   void (*int_derivs)(), int npts, 
   double *z, double *fnz, double *fnzspline);
extern void myrkqc(
   double *y, double *dydx, double *x, double htry, double eps,
   double *yscal, double *hdid, double *hnext,
   void (*int_derivs)(), int npts, 
   double *z, double *fnz, double *fnzspline);
extern void rknstp(
   double *y, double *yp, double *d2ydx2, int nvar, 
   double *x, double h, double eps,
   double *hdid, double *hnext,
   void (*derivs)(), double *odeparams, int *flags,
   double zind[], double U[], double U2[], double Upp[], double Upp2[],
   double nsqr[], double nsqr2[], int nind);
extern void rknde(
   double *ystart, double *ypstart, int nvar, double xin, double xout, 
   double eps, int *nok, int *nbad,
   void (*derivs)(), double *odeparams, int *flags,
   double zind[], double U[], double U2[], double Upp[], double Upp2[],
   double nsqr[], double nsqr2[], int nind,
   int printflg, double *normphi, double **efn, int *tnstp);
extern void myerny(
   double *phi, int nvar, double xin, double xout, 
   double eps, int *nok, int *nbad,
   void (*derivs)(), double *odeparams, int *flags,
   double zind[], double U[], double U2[], double Upp[], double Upp2[],
   double nsqr[], double nsqr2[], int nind,
   int printflg, double *normphi, double **efn, int *tnstp);

extern double integrate(
   double *fx, double *x, int npts, double x1, double x2);
extern double intfdx(
   double *x, double *f, int npts, double x1, double x2);
extern void intf(
   double *f, double *x, int npts, double *ifn);
extern void integrateO0(
   double *f, int ox, int nx, double dx, double if0, double *ifn);
extern void integrateO1(
   double *f, int ox, int nx, double dx, int isgn, double if0, double *ifn);
extern void integrateO1_gen(
   double *f, double *x, int ox, int nx, int isgn, double if0, double *ifn);

extern void CintegrateO0(
   dcomplex *f, int ox, int nx, double dx, dcomplex if0, dcomplex *ifn);
extern void CintegrateO0_gen(
   dcomplex *f, double *x, int ox, int nx, int isgn, 
   dcomplex if0, dcomplex *ifn);
extern void CintegrateO1(
   dcomplex *f, int ox, int nx, double dx, int isgn, 
   dcomplex if0, dcomplex *ifn);
extern void CintegrateO1_gen(
   dcomplex *f, double *x, int ox, int nx, int isgn, 
   dcomplex if0, dcomplex *ifn);



extern void ZintegrateO0(
   zomplex *f, int ox, int nx, double dx, zomplex if0, zomplex *ifn);
extern void ZintegrateO0_gen(
   zomplex *f, double *x, int ox, int nx, int isgn, 
   zomplex if0, zomplex *ifn);
extern void ZintegrateO1(
   zomplex *f, int ox, int nx, double dx, int isgn, 
   zomplex if0, zomplex *ifn);
extern void ZintegrateO1_gen(
   zomplex *f, double *x, int ox, int nx, int isgn, 
   zomplex if0, zomplex *ifn);
