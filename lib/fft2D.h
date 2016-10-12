extern void hfft2D(
   double **mat, int nn, int mm, dcomplex **cfld, int isgn);
extern void vrfft2D(
   double **mat, int nn, int mm, dcomplex **cfld, int isgn);
extern void vfft2D(
   dcomplex **cf1, int nn, int mm, dcomplex **cf2, int isgn);
extern void vsfft2D(
   dcomplex **cf1, int nn, int mm, dcomplex **cf2, int isgn);
extern void vcfft2D(
   dcomplex **cf1, int nn, int mm, dcomplex **cf2, int isgn);
extern void sfft2D(
   double **m1, int nn, int mm, int index, double **m2, int isgn);
extern void cfft2D(
   double **m1, int nn, int mm, int index, double **m2, int isgn);
extern void znfft2D(
   double **mat, int nn, int om, int mm, dcomplex **cfld, int isgn);
extern void azfft2D(
   double **mat, int on, int nn, int mm, dcomplex **cfld, int isgn);

extern void zhfft2D(
   double **mat, int nn, int mm, zomplex **cfld, int isgn);
extern void zvrfft2D(
   double **mat, int nn, int mm, zomplex **cfld, int isgn);
extern void zvfft2D(
   zomplex **cf1, int nn, int mm, zomplex **cf2, int isgn);
extern void zvsfft2D(
   zomplex **cf1, int nn, int mm, zomplex **cf2, int isgn);
extern void zvcfft2D(
   zomplex **cf1, int nn, int mm, zomplex **cf2, int isgn);
extern void zznfft2D(
   double **mat, int nn, int om, int mm, zomplex **cfld, int isgn);
extern void zazfft2D(
   double **mat, int on, int nn, int mm, zomplex **cfld, int isgn);

