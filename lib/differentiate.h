extern double intplderiv1(
   double x1, double x2,
   double y1, double y2);

extern double intplderiv(
   double xv, 
   double x1, double x2, double x3,
   double y1, double y2, double y3);

extern void diff(double *fn, double *x, int nind, double *fnp);
extern void diffO1(double *fn, double *x, int nind, double *fnp, int isign);
extern void diff1D(
   double *f, int ox, int nx, double dx, double *df);
extern void Pdiff1D(
   double *f, int ox, int nx, double dx, double *df);

extern void pdiff2DO1(
   double **f, int ox, int nx, int oy, int ny,
   double dxy, int index, double **df);
extern void Ppdiff2DO1(
   double **f, int ox, int nx, int oy, int ny,
   double dxy, int index, double **df);
extern void Zpdiff2DO1(
   double **f, int ox, int nx, int oy, int ny,
   double dxy, int index, double **df);


extern void pdiff2D(
   double **f, int ox, int nx, int oy, int ny,
   double dxy, int index, double **df);
extern void pdiff2D_gen(
   double **f, double *x, int ox, int nx, double *y, int oy, int ny,
   int index, double **df);
extern void Ppdiff2D(
   double **f, int ox, int nx, int oy, int ny,
   double dxy, int index, double **df);


extern void Cpdiff2DO1(
   dcomplex **f, int ox, int nx, int oy, int ny,
   double dxy, int index, dcomplex **df);
extern void CPpdiff2DO1(
   dcomplex **f, int ox, int nx, int oy, int ny,
   double dxy, int index, dcomplex **df);
extern void Cpdiff2D(
   dcomplex **f, int ox, int nx, int oy, int ny,
   double dxy, int index, dcomplex **df);
extern void Cpdiff2D_gen(
   dcomplex **f, 
   double *x, int ox, int nx, 
   double *y, int oy, int ny,
   int index, dcomplex **df);
extern void CPpdiff2D(
   dcomplex **f, int ox, int nx, int oy, int ny,
   double dxy, int index, dcomplex **df);
extern void CFpdiff2D(
   dcomplex **f, int ox, int nx, int oy, int ny,
   double *kxy, int index, dcomplex **df);
extern void CFpdiff2D_gen(
   dcomplex **f, int ox, int nx, int oy, int ny,
   double *kxy, int index, dcomplex **df);

extern void Zzpdiff2DO1(
   zomplex **f, int ox, int nx, int oy, int ny,
   double dxy, int index, zomplex **df);
extern void ZPpdiff2DO1(
   zomplex **f, int ox, int nx, int oy, int ny,
   double dxy, int index, zomplex **df);
extern void Zpdiff2D(
   zomplex **f, int ox, int nx, int oy, int ny,
   double dxy, int index, zomplex **df);
extern void Zpdiff2D_gen(
   zomplex **f, 
   double *x, int ox, int nx, 
   double *y, int oy, int ny,
   int index, zomplex **df);
extern void ZPpdiff2D(
   zomplex **f, int ox, int nx, int oy, int ny,
   double dxy, int index, zomplex **df);
extern void ZFpdiff2D(
   zomplex **f, int ox, int nx, int oy, int ny,
   double *kxy, int index, zomplex **df);
extern void ZFpdiff2D_gen(
   zomplex **f, int ox, int nx, int oy, int ny,
   double *kxy, int index, zomplex **df);

extern void pdiff3DO1(
   double ***f, int ox, int nx, int oy, int ny, int oz, int nz,
   double dxyz, int index, double ***df);
extern void Ppdiff3DO1(
   double ***f, int ox, int nx, int oy, int ny, int oz, int nz,
   double dxyz, int index, double ***df);
extern void Zpdiff3DO1(
   double ***f, int ox, int nx, int oy, int ny, int oz, int nz,
   double dxyz, int index, double ***df);

extern void Cpdiff3DO1(
   dcomplex ***f, int ox, int nx, int oy, int ny, int oz, int nz,
   double dxyz, int index, dcomplex ***df);
extern void CPpdiff3DO1(
   dcomplex ***f, int ox, int nx, int oy, int ny, int oz, int nz,
   double dxyz, int index, dcomplex ***df);
extern void Cpdiff3D(
   dcomplex ***f, int ox, int nx, int oy, int ny, int oz, int nz,
   double dxyz, int index, dcomplex ***df);
extern void Cpdiff3D_gen(
   dcomplex ***f, 
   double *x, int ox, int nx, 
   double *y, int oy, int ny, 
   double *z, int oz, int nz,
   int index, dcomplex ***df);
extern void CFpdiff3D(
   dcomplex ***f, int ox, int nx, int oy, int ny, int oz, int nz,
   double *kxyz, int index, dcomplex ***df);


extern void Zzpdiff3DO1(
   zomplex ***f, int ox, int nx, int oy, int ny, int oz, int nz,
   double dxyz, int index, zomplex ***df);
extern void ZPpdiff3DO1(
   zomplex ***f, int ox, int nx, int oy, int ny, int oz, int nz,
   double dxyz, int index, zomplex ***df);
extern void Zpdiff3D(
   zomplex ***f, int ox, int nx, int oy, int ny, int oz, int nz,
   double dxyz, int index, zomplex ***df);
extern void Zpdiff3D_gen(
   zomplex ***f, 
   double *x, int ox, int nx, 
   double *y, int oy, int ny, 
   double *z, int oz, int nz,
   int index, zomplex ***df);
extern void ZFpdiff3D(
   zomplex ***f, int ox, int nx, int oy, int ny, int oz, int nz,
   double *kxyz, int index, zomplex ***df);
