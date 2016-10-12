#include <math.h>
#include "complex.h"
#include "mylaguer.h"
#include "constants.h"
#include "macros.h"

#define EPS 2.0e-6
#define MAXM 100

void myzroots(a,m,roots,polish)
   int m,polish;
   dcomplex *a,*roots;
{
   int jj,j,i;
   dcomplex x,b,c,ad[MAXM];
   dcomplex cZero;

   cZero=Complex(0.0,0.0);

   for (j=0;j<=m;j++) ad[j]=a[j];
   for (j=m;j>=1;j--) {
      x=cZero;
      mylaguer(ad,j,&x,EPS,0);
      if (dabs(x.i) <= (2.0*EPS*dabs(x.r))) x.i=0.0;
      roots[j]=x;
      b=ad[j];
      for (jj=j-1;jj>=0;jj--) {
         c=ad[jj];
         ad[jj]=b;
         b=Cadd(Cmul(x,b),c);
      }
   }
   if (polish)
      for (j=1;j<=m;j++)
         mylaguer(a,m,&roots[j],EPS,1);
   for (j=2;j<=m;j++) {
      x=roots[j];
      for (i=j-1;i>=1;i--) {
         if (roots[i].r <= x.r) break;
         roots[i+1]=roots[i];
      }
      roots[i+1]=x;
   }
}

#undef EPS
#undef MAXM
