extern void scaleshortfldO1(
   short **f1, int n1, int m1,
   short **f2, int n2,int m2);
extern void scalefldO1(
   double **f1, int n1, int m1,
   double **f2, int n2, int m2);
extern void scalefn(double *x1,double *f1,int n1, 
                    double *x2,double *f2,int n2);
extern void scalefld(
   double **f1, int n1, int m1, int i0flg,
   double **f2, int n2, int m2, int o0flg);
extern void scalefield(
   double **f1, int ox1, int nx1, int oy1, int ny1, 
   double xmn1, double xmx1, double ymn1, double ymx1,
   double **f2, int ox2, int nx2, int oy2, int ny2, 
   double xmn2, double xmx2, double ymn2, double ymx2,
   int exflg);
extern void rotatefldO1(
   double **f1, int ox1, int nx1, int oy1, int ny1, 
   double xmn1, double xmx1, double ymn1, double ymx1,
   double Ox, double Oy, double rot,
   double **f2,
   int exflg);
extern void reflectfldO1(
   double **f1, int ox1, int nx1, int oy1, int ny1, 
   double xmn1, double xmx1, double ymn1, double ymx1,
   double x1, double y1, double x2, double y2,
   double **f2,
   int exflg);
extern void copyfield(
   double **f1, int ox, int nx, int oy, int ny, 
   double **f2);
