extern void Rstagger(
   double *f, int ox, int nx, int staggerflg, double *sf
);
extern void _Rstagger(
   double *f, int ox, int nx, int staggerflg
);
extern void Cstagger(
   dcomplex *f, int ox, int nx, int staggerflg, dcomplex *sf
);

extern void Rstagger2D(
   double **f, int oy, int ny, int oz, int nz, int staggerflg, double **sf
);
extern void _Rstagger2D(
   double **f, int oy, int ny, int oz, int nz, int staggerflg
);
extern void Cstagger2D(
   dcomplex **f, int oy,int ny, int oz,int nz, int staggerflg, dcomplex **sf
);
extern void _Cstagger2D(
   dcomplex **f, int oy,int ny, int oz,int nz, int staggerflg
);

extern void Rstagger3D(
   double ***f, int ox,int nx, int oy,int ny, int oz,int nz, 
   int staggerflg, double ***sf
);
extern void _Rstagger3D(
   double ***f, int ox,int nx, int oy,int ny, int oz,int nz, 
   int staggerflg
);
extern void Cstagger3D(
   dcomplex ***f, int ox,int nx, int oy,int ny, int oz,int nz, 
   int staggerflg, dcomplex ***sf
);
extern void _Cstagger3D(
   dcomplex ***f, int ox,int nx, int oy,int ny, int oz,int nz, 
   int staggerflg
);
extern void RPstagger2D(
   double **f, int ox,int nx, int oy,int ny,
   int staggerflg, double **sf
);
extern void _RPstagger2D(
   double **f, int ox,int nx, int oy,int ny,
   int staggerflg
);
extern void RPstagger3D(
   double ***f, int ox,int nx, int oy,int ny, int oz,int nz, 
   int staggerflg, double ***sf
);
extern void _RPstagger3D(
   double ***f, int ox,int nx, int oy,int ny, int oz,int nz, 
   int staggerflg
);
extern void RZstagger2D(
   double **f, int ox,int nx, int oy,int ny,
   int staggerflg, double **sf
);
extern void _RZstagger2D(
   double **f, int ox,int nx, int oy,int ny,
   int staggerflg
);
extern void RZstagger3D(
   double ***f, int ox,int nx, int oy,int ny, int oz,int nz, 
   int staggerflg, double ***sf
);
extern void _RZstagger3D(
   double ***f, int ox,int nx, int oy,int ny, int oz,int nz, 
   int staggerflg
);
