extern double  myzbrent(
   double (*func)(double x), double x1, double x2, double tol);

extern double  myzbrent2(
   double (*func)(int n, double x), int n, double x1, double x2, double tol);
