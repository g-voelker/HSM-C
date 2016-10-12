extern void _Ccmul2D(
   dcomplex **f, double c, int ox, int nx, int oy, int ny);

extern void _Cadd2D(
   dcomplex **f1, dcomplex **f2, int ox, int nx, int oy, int ny);

extern void _Csub2D(
   dcomplex **f1, dcomplex **f2, int ox, int nx, int oy, int ny);

extern void Radd3D(
   double ***f1, double ***f2, 
   int ox, int nx, int oy, int ny, int oz, int nz,
   double ***f3);

extern void _Radd3D(
   double ***f1, double ***f2, 
   int ox, int nx, int oy, int ny, int oz, int nz
);

extern void Rsub3D(
   double ***f1, double ***f2, 
   int ox, int nx, int oy, int ny, int oz, int nz,
   double ***f3);

extern void _Rsub3D(
   double ***f1, double ***f2, 
   int ox, int nx, int oy, int ny, int oz, int nz
);

extern void Rmul3D(
   double ***f1, double ***f2, 
   int ox, int nx, int oy, int ny, int oz, int nz,
   double ***f3);

extern void _Rmul3D(
   double ***f1, double ***f2, 
   int ox, int nx, int oy, int ny, int oz, int nz
);

extern void Rdiv3D(
   double ***f1, double ***f2,
   int ox, int nx, int oy, int ny, int oz, int nz,
   double ***f3);

extern void _Rdiv3D(
   double ***f1, double ***f2,
   int ox, int nx, int oy, int ny, int oz, int nz);

extern void Rcmul3D(
   double ***f1, double c, int ox, int nx, int oy, int ny, int oz, int nz,
   double ***f2);

extern void _Rcmul3D(
   double ***f1, double c, int ox, int nx, int oy, int ny, int oz, int nz);

extern void Rmul3D1D(
   double ***f1, double *f2, int index,
   int ox, int nx, int oy, int ny, int oz, int nz,
   double ***f3);

extern void _Rmul3D1D(
   double ***f1, double *f2, int index,
   int ox, int nx, int oy, int ny, int oz, int nz);

extern void Rmul3D2D(
   double ***f1, double **f2, int index,
   int ox, int nx, int oy, int ny, int oz, int nz,
   double ***f3);

extern void _Rmul3D2D(
   double ***f1, double **f2, int index,
   int ox, int nx, int oy, int ny, int oz, int nz);

extern void Rdiv3D2D(
   double ***f1, double **f2, int index,
   int ox, int nx, int oy, int ny, int oz, int nz,
   double ***f3);

extern void _Rdiv3D2D(
   double ***f1, double **f2, int index,
   int ox, int nx, int oy, int ny, int oz, int nz);

extern void _Ccmul3D(
   dcomplex ***f, double c, 
   int ox, int nx, int oy, int ny, int oz, int nz);

extern void _CRvmul3D(
   dcomplex ***f, double *v, int index,
   int ox, int nx, int oy, int ny, int oz, int nz);

extern void _Cadd3D(
   dcomplex ***f1, dcomplex ***f2,
   int ox, int nx, int oy, int ny, int oz, int nz);

extern void _Csub3D(
   dcomplex ***f1, dcomplex ***f2,
   int ox, int nx, int oy, int ny, int oz, int nz);

extern void Rmbsmooth2D(
   double **f, int ox, int nx, int oy, int ny, double **sf);
extern void _Rmbsmooth2D(
   double **f, int ox, int nx, int oy, int ny);
extern void Rmbsmooth(
   double ***f, int ox, int nx, int oy, int ny, int os, int ns, double ***sf);
extern void _Rmbsmooth(
   double ***f, int ox, int nx, int oy, int ny, int os, int ns);

extern void _Cmbsmooth2D(
   dcomplex **f, int ox, int nnx, int oy, int ny);
extern void _Cmbsmooth2D_2(
   dcomplex **f, int ox, int nnx, int oy, int ny);
extern void _Cmbsmooth(
   dcomplex ***f, int ox, int nnx, int oy, int ny, int os, int ns);

extern void _Rmbzero2D(
   double **f, int ox, int nnx, int oy, int ny);

extern void _Rmbzero(
   double ***f, int ox, int nnx, int oy, int ny, int os, int ns);

extern void _Cmbzero2D(
   dcomplex **f, int ox, int nnx, int oy, int ny);

extern void _Cmbzero(
   dcomplex ***f, int ox, int nnx, int oy, int ny, int os, int ns);
