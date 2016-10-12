#include <stdio.h>
#include <math.h>
#include "macros.h"

#define ITMAX 100    /* maximum number of iterations */
#define EPS 5.0e-16  /* machine double precision */

/* =============================================================================
   MYZBRENT: Sept 28/1995

   Zero finder, adapted from numerical recipes.
   =============================================================================
*/
double myzbrent(double (*func)(), double x1, double x2, double tol)
{
   int iter;
   double a=x1,b=x2,c,d,e,min1,min2;
   double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s, tol1,xm; 

   fa=(*func)(a);
   fb=(*func)(b);

   if(fb*fa>0.0) {
       printf("MYZBRENT: (a,fa)=(%10.3e,%10.3e) (b,fb)=(%10.3e,%10.3e)\n",
          a,fa,b,fb); 
       myerror("MYZBRENT (ERROR): Root must be bracketed");
   }
   fc=fb;
   for (iter=1;iter<=ITMAX;iter++) {
      if(fb*fc>0.0) {
         c=a;
         fc=fa;
         e=d=b-a;
      }
      if(dabs(fc)<dabs(fb)) {
         a=b;
         b=c;
         c=a;
         fa=fb;
         fb=fc;
         fc=fa;
      }
      tol1=2.*EPS*dabs(b)+0.5*tol;
      xm=0.5*(c-b);
      if(dabs(xm)<=tol1 || fb==0.0) return b;
      if(dabs(e)>=tol1 && dabs(fa)>dabs(fb)) {
         s=fb/fa;
         if(a==c) {
            p=2.*xm*s;
            q=1.0-s;
         } else {
            q=fa/fc;
            r=fb/fc;
            p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
            q=(q-1.0)*(r-1.0)*(s-1.0);
         }
         if(p>0.0) q=-q;
         p=dabs(p);
         min1 = 3.0*xm*q-dabs(tol1*q);
         min2 = dabs(e*q);
         if(2.0*p < (min1<min2 ? min1 : min2)) {
            e=d;
            d=p/q;
         } else {
            d=xm;
            e=d;
         }
      } else {
         d=xm;
         e=d;
      }
      a=b;
      fa=fb;
      if(dabs(d) > tol1) b += d;
      else b += (xm > 0.0 ? dabs(tol1) : -dabs(tol1));
      fb= (*func)(b);
   }
   myerror("ZBRENT (ERROR): exceeding maximum iterations");
   return 0.0;
}

/* =============================================================================
   MYZBRENT2: June 9/2012

   Zero finder for function with integer index as second variables
   (Needed for zeros of nth order bessel function).

   =============================================================================
*/
double myzbrent2(double (*func)(int n, double x), 
                 int n, double x1, double x2, double tol)
{
   int iter;
   double a=x1,b=x2,c,d,e,min1,min2;
   double fa,fb,fc,p,q,r,s, tol1,xm; 

   fa=(*func)(n,a);
   fb=(*func)(n,b);

   if(fb*fa>0.0) {
       printf("MYZBRENT2: (a,fa)=(%10.3e,%10.3e) (b,fb)=(%10.3e,%10.3e)\n",
          a,fa,b,fb); 
       myerror("MYZBRENT2 (ERROR): Root must be bracketed");
   }
   fc=fb;
   for (iter=1;iter<=ITMAX;iter++) {
      if(fb*fc>0.0) {
         c=a;
         fc=fa;
         e=d=b-a;
      }
      if(dabs(fc)<dabs(fb)) {
         a=b;
         b=c;
         c=a;
         fa=fb;
         fb=fc;
         fc=fa;
      }
      tol1=2.*EPS*dabs(b)+0.5*tol;
      xm=0.5*(c-b);
      if(dabs(xm)<=tol1 || fb==0.0) return b;
      if(dabs(e)>=tol1 && dabs(fa)>dabs(fb)) {
         s=fb/fa;
         if(a==c) {
            p=2.*xm*s;
            q=1.0-s;
         } else {
            q=fa/fc;
            r=fb/fc;
            p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
            q=(q-1.0)*(r-1.0)*(s-1.0);
         }
         if(p>0.0) q=-q;
         p=dabs(p);
         min1 = 3.0*xm*q-dabs(tol1*q);
         min2 = dabs(e*q);
         if(2.0*p < (min1<min2 ? min1 : min2)) {
            e=d;
            d=p/q;
         } else {
            d=xm;
            e=d;
         }
      } else {
         d=xm;
         e=d;
      }
      a=b;
      fa=fb;
      if(dabs(d) > tol1) b += d;
      else b += (xm > 0.0 ? dabs(tol1) : -dabs(tol1));
      fb= (*func)(n,b);
   }
   myerror("MYZBRENT2 (ERROR): exceeding maximum iterations");
   return 0.0;
}
