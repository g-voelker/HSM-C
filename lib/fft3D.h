extern void znfft3D(
   double ***f, 
   int nx, int oy, int ny, int oz, int nz, dcomplex ***F, int isgn);
extern void mdfft3D(
   dcomplex ***f, 
   int nx, int ny, int oz, int nz, dcomplex ***F, int isgn);
extern void vfft3D(
   dcomplex ***f, 
   int nx, int ny, int nz, dcomplex ***F, int isgn);

extern void zznfft3D(
   double ***f, 
   int nx, int oy, int ny, int oz, int nz, zomplex ***F, int isgn);
extern void zmdfft3D(
   zomplex ***f, 
   int nx, int ny, int oz, int nz, zomplex ***F, int isgn);
extern void zvfft3D(
   zomplex ***f, 
   int nx, int ny, int nz, zomplex ***F, int isgn);
