/* ============================================================================
   Files to perform statistics and simple analyses on functions and fields
   ============================================================================
*/
#include "complex.h"
#include "alloc_space.h"
#include "constants.h"
#include "macros.h"

/* =========================================================================
   Find extremum and its location of function evaluated at 3 points 
   (x1,y1), (x2,y2), (x3,y3)
   =========================================================================
*/
int getextremum(
  double x1, double x2, double x3, double y1, double y2, double y3,
  double *xm, double *ym)
{
  double m21,m32,A,B,C;

  if (x1==x2 || x1==x3 || x2==x3) return(-1);

  m21 = (y2-y1)/(x2-x1);
  m32 = (y3-y2)/(x3-x2);
  A = (m21-m32)/(x1-x3);
  B = m21 - (x1+x2)*A;
  C = y1 - (A*x1 + B)*x1;

  if (A==0) return(-2);

  *xm = -0.5*B/A;
  *ym = C-0.25*B*B/A;

  if (A>0) { /* minimum */
    return(0);
  } else {   /* maximum */
    return(1);
  }
}

/* =====================================================================
   Find autocorrelation of function evaluated at a particular index.

   Autocorrelation is given by: R(m) = sum_m^N (f[n] f[n-m])
   =====================================================================
*/
double indexautocorrel(
  int m, double *f, int N)
{
  int n;
  double cf;

  cf = 0.0;
  for(n=0;n<=N;n++) cf += f[n]*f[n-m];

  return(cf);
}

/* =====================================================================
   Find correlation of two functions evaluated at a particular index.

   Correlation of two functions is given by: R(m) = sum_n (f1[n] f2[n-m])

   if the indices of f1 and f2 range from 0..N, then

                   R(m) = sum_(n=nmin..nmax)   f1[n] f2[n-m]

   where nmin = max(0,m) and nmax=min(N,N+m).
   =====================================================================
*/
double indexcorrel(
  int m, double *f0, double *f1, int N)
{
  int n,nmin,nmax;
  double cf;

  nmin=max(0,m);
  nmax=min(N,N+m);

  cf = 0.0;
  for(n=nmin;n<=nmax;n++) cf += f0[n]*f1[n-m];

  return(cf);
}
