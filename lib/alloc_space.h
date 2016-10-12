extern void myerror(char *);
extern short *svector(int o1, int n1);
extern void free_svector(short *iv, int o1, int n1);
extern int *ivector(int o1, int n1);
extern void free_ivector(int *iv, int o1, int n1);
extern float *fvector(int o1, int n1);
extern void free_fvector(float *fv,int o1, int n1);
extern double *dvector(int o1, int n1);
extern void free_dvector(double *dv,int o1, int n1);
extern dcomplex *dcvector(int o1, int n1);
extern void free_dcvector(dcomplex *dcv,int o1, int n1);
extern zomplex *zvector(int o1, int n1);
extern void free_zvector(zomplex *dcv,int o1, int n1);
extern short **smatrix2(int o1, int n1, int o2, int n2);
extern void free_smatrix2(short **im,int o1, int n1, int o2, int n2);
extern int **imatrix2(int o1, int n1, int o2, int n2);
extern void free_imatrix2(int **im,int o1, int n1, int o2, int n2);
extern float **fmatrix2(int o1, int n1, int o2, int n2);
extern void free_fmatrix2(float **fm,int o1, int n1, int o2, int n2);
extern double **dmatrix2(int o1, int n1, int o2, int n2);
extern void free_dmatrix2(double **dm,int o1, int n1, int o2, int n2);
extern dcomplex **dcmatrix2(int o1, int n1, int o2, int n2);
extern void free_dcmatrix2(dcomplex **dcm,int o1, int n1, int o2, int n2);
extern zomplex **zmatrix2(int o1, int n1, int o2, int n2);
extern void free_zmatrix2(zomplex **dcm,int o1, int n1, int o2, int n2);
extern int ***imatrix3(int o1, int n1, int o2, int n2,int o3, int n3);
extern void free_imatrix3(
   int ***im,
   int o1, int n1, int o2, int n2,int o3, int n3);
extern short ***smatrix3(int o1, int n1, int o2, int n2,int o3, int n3);
extern void free_smatrix3(
   short ***im,
   int o1, int n1, int o2, int n2,int o3, int n3);
extern float ***fmatrix3(
   int o1, int n1, int o2, int n2,int o3, int n3);
extern void free_fmatrix3(
   float ***fm,
   int o1, int n1, int o2, int n2,int o3, int n3);
extern double ***dmatrix3(
   int o1, int n1, int o2, int n2,int o3, int n3);
extern void free_dmatrix3(
   double ***dm, 
   int o1, int n1, int o2, int n2,int o3, int n3);
extern dcomplex ***dcmatrix3(
   int o1, int n1, int o2, int n2,int o3, int n3);
extern void free_dcmatrix3(
   dcomplex ***dcm,
   int o1, int n1, int o2, int n2,int o3, int n3);
extern zomplex ***zmatrix3(
   int o1, int n1, int o2, int n2,int o3, int n3);
extern void free_zmatrix3(
   zomplex ***dcm,
   int o1, int n1, int o2, int n2,int o3, int n3);
extern float ****fmatrix4(
   int o1, int n1, int o2, int n2,int o3, int n3,int o4, int n4);
extern void free_fmatrix4(
   float ****fm,
   int o1, int n1, int o2, int n2,int o3, int n3,int o4, int n4);
extern double ****dmatrix4(
   int o1, int n1, int o2, int n2,int o3, int n3,int o4, int n4);
extern void free_dmatrix4(
   double ****dm,
   int o1, int n1, int o2, int n2,int o3, int n3,int o4, int n4);
extern dcomplex ****dcmatrix4( 
   int o1, int n1, int o2, int n2,int o3, int n3,int o4, int n4);
extern void free_dcmatrix4(
   dcomplex ****dcm,
   int o1, int n1, int o2, int n2,int o3, int n3,int o4, int n4);
extern zomplex ****zmatrix4( 
   int o1, int n1, int o2, int n2,int o3, int n3,int o4, int n4);
extern void free_zmatrix4(
   zomplex ****dcm,
   int o1, int n1, int o2, int n2,int o3, int n3,int o4, int n4);
