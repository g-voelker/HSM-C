/* Hyperbolic Functions */
extern double sech(double x);
extern double sechsqr(double x);
extern double tanh(double x);
extern double tanhsqr(double x);

/* Bessel Functions */

extern double mybessj0(double r);
extern double mybessj1(double r);
extern double mybessj2(double r);
extern double mybessjn(int n,double r);
extern double mybessjnhalf(int n,double r);
extern double mybessy0(double r);
extern double mybessy1(double r);
extern double mybessyn(int n,double r);

extern double mydbessj0(double r);
extern double mydbessj1(double r);
extern double mydbessjn(int n,double r);

extern double bessj0zero(int m);
extern double bessj1zero(int m);
extern void bessjnzeros(int n, int nzeros, double *Jz);

extern double mybessi0(double r);
extern double mybessi1(double r);
extern double mybessk0(double r);
extern double mybessk1(double r);
extern double mybessin(int n,double r);
extern double mybesskn(int n,double r);

/* Legendre Polynomials */
extern double pn(int n,double x);
extern double pnm(int n,int m,double x);

/* Hermite Polynomials */
extern double hn(int n,double x);

/* Frenel Sine and Cosine Integrals */
extern void frenel(double x, double *s, double *c);

/* Dawson Integrals */
extern double dawson(double x);

