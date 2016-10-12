extern void defj0array(
  double **J0array, int nr, int nz, double *J0z);

extern void nrmJ0ser(
  double *nrmJ0,int nz,double *J0z);

extern void j0trans(
  double *fnr,double *fnbc,double **J0array, int nr,int nz,
  double *J0z,double *nrmJ0);

extern void j0invtrans(
  double *fnr,double *fnbc,double **J0array,int nr,int nz);

extern void j0cfseries(
  double *fnr,int nr, double *fnbc,int nz);

extern void j0sumseries(
  double *fnbc,int nz,double *fnr,int nr);

extern void defj1array(
  double **J1array, int nr, int nz, double *J1z);

extern void nrmJ1ser(
  double *nrmJ1,int nz,double *J1z);

extern void j1trans(
  double *fnr,double *fnbc,double **J1array, int nr,int nz,
  double *J1z,double *nrmJ1);

extern void j1invtrans(
  double *fnr,double *fnbc,double **J1array,int nr,int nz);

extern void j1cfseries(
  double *fnr,int nr, double *fnbc,int nz);

extern void j1sumseries(
  double *fnbc,int nz,double *fnr,int nr);


extern void defjnarray(
  double **JNarray, int n, int nr, int nz, double *JNz);

extern void nrmJNser(
  double *nrmJN,int n,int nz,double *JNz);

extern void jtrans(
  double *fnr,double *fnbc,double **JNarray, int nr,int nz,
  double *JNz,double *nrmJN);

extern void jinvtrans(
  double *fnr,double *fnbc,double **JNarray,int nr,int nz);

