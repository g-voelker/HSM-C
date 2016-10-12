extern void myfour1(
   double *data,int nn,int isign);

extern void myrealft(
   double *data,int n,int isign);

extern void mycosft(
   double *y,int n);

extern void mycosft2(
   double *y, int n, int isign);

extern void mysinft(
   double *y,int n);

extern void mytwofft(
   double *data1,double *data2,double *fft1,double *fft2,int n);

extern void myconvlv(
   double *data,int n,double *respns,int m,int rfflg,
   int isign,double *ans);

extern void mycorrel(
   double *data1,double *data2,int n, double *ans);

extern void myfft(
   double *data,int nn, dcomplex *fdata, int isign);

extern void myzfft(
   double *data,int nn, zomplex *fdata, int isign);

extern void mycfft(
   dcomplex *data,int mm, dcomplex *fdata, int isign);

extern void myczfft(
   zomplex *data, int mm, zomplex *fdata, int isign);
