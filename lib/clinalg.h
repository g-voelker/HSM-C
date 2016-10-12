typedef dcomplex **cmatrix_t, *cvector_t;
typedef zomplex **zmatrix_t, *zvector_t;

extern void ludcmp(
   dcomplex **a, int n, int *indx, dcomplex *d);
extern void lubksb(
   dcomplex **a, int n, int *indx, dcomplex d[]);
extern void printCmat(cmatrix_t c, int n);
extern void Cmatvectmult(cmatrix_t a, cvector_t b, cvector_t c, int n);
extern void Cmatmult(cmatrix_t a, cmatrix_t b, cmatrix_t c, int n);
extern dcomplex Cdotprod(cvector_t a,cvector_t b, int n);
extern void Cinvert(cmatrix_t a,cmatrix_t b, int n);

extern void zludcmp(
   zomplex **a, int n, int *indx, zomplex *d);
extern void zlubksb(
   zomplex **a, int n, int *indx, zomplex d[]);
extern void printZmat(zmatrix_t c, int n);
extern void Zmatvectmult(zmatrix_t a, zvector_t b, zvector_t c, int n);
extern void Zmatmult(zmatrix_t a, zmatrix_t b, zmatrix_t c, int n);
extern zomplex Zdotprod(zvector_t a,zvector_t b, int n);
extern void Zinvert(zmatrix_t a,zmatrix_t b, int n);
