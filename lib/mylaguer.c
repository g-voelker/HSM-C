#include <math.h>
#include "myerror.h"
#include "complex.h"
#include "constants.h"
#include "macros.h"

#define EPSS 6.e-8
#define MAXIT 500

void mylaguer(a,m,x,eps,polish)
   int m,polish;
   double eps;
   dcomplex *a,*x;
{
   int j,iter;
   double err,dxold,cdx,abx;
   dcomplex sq,h,gp,gm,g2,g,b,d,dx,f,x1;
   dcomplex cZero;

   cZero = Complex(0.0,0.0);
   dxold=Cabs(*x);
   for (iter=1;iter<=MAXIT;iter++) {
      b=a[m];
      err=Cabs(b);
      d=f=cZero;
      abx=Cabs(*x);
      for (j=m-1;j>=0;j--) {
         f=Cadd(Cmul(*x,f),d);
         d=Cadd(Cmul(*x,d),b);
         b=Cadd(Cmul(*x,b),a[j]);
         err=Cabs(b)+abx*err;
      }
      err *= EPSS;
      if (Cabs(b) <= err) return;
      g=Cdiv(d,b);
      g2=Cmul(g,g);
      h=Csub(g2,RCmul(2.0,Cdiv(f,b)));
      sq=Csqrt(RCmul((double) (m-1),Csub(RCmul((double) m,h),g2)));
      gp=Cadd(g,sq);
      gm=Csub(g,sq);
      if (Cabs(gp) < Cabs(gm))gp=gm;
      dx=Cdiv(Complex((double) m,0.0),gp);
      x1=Csub(*x,dx);
      if (x->r == x1.r && x->i == x1.i) return;
      *x=x1;
      cdx=Cabs(dx);
      if (iter > 6 && cdx >= dxold) return;
      dxold=cdx;
      if (!polish)
         if (cdx <= eps*Cabs(*x)) return;
   }
   myerror("MYLAGUER: Too many iterations");
}

#undef EPSS
#undef MAXIT
