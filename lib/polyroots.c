/* ==================================================================
   POLYROOTS: Nov. 27, 2003
   Find roots of polynomials given coefficients 
   ==================================================================
*/
#include <math.h>
#include "complex.h"
#include "alloc_space.h"
#include "constants.h"
#include "myerror.h"
#include "macros.h"

#define MINREAL 1.0e-10

void solvecubic(double *cf, dcomplex *root)
{
   int i;
   double p,q,D,sqrtD,y1cubed,y2cubed,y1,y2,A3,theta3,A,theta;
   double scf[4];
   dcomplex c3rtone,cycubed,cy1,cy2;

   c3rtone = Complex(-0.5,0.5*Root3);
 
   for (i=0;i<=3;i++) scf[i] = cf[i]/cf[3];

   /* make substitution x -> y-b/3 */
   p = -scf[2]*scf[2]/3.0 + scf[1];
   q = 2*scf[2]*scf[2]*scf[2]/27.0 -scf[2]*scf[1]/3.0 + scf[0];

   D = q*q/4.0 + p*p*p/27.0;
   if (D>=0.0) {
      sqrtD = sqrt(D);
      y1cubed = -q/2.0 + sqrtD;
      y2cubed = -q/2.0 - sqrtD;

      y1cubed = sgn(y1cubed)*max(MINREAL,dabs(y1cubed));
      y2cubed = sgn(y2cubed)*max(MINREAL,dabs(y2cubed));
      y1 = sgn(y1cubed)*exp(log(dabs(y1cubed))/3.0);
      y2 = sgn(y2cubed)*exp(log(dabs(y2cubed))/3.0);

      root[1].r = y1 + y2 - scf[2]/3.0;  root[1].i=0.0;
      root[2].r = -(y1 + y2)/2.0 - scf[2]/3.0;  root[2].i = Root3*(y1 - y2)/2.0;
      root[3] = root[2];  root[3].i *= -1.0;
   } else {
      sqrtD = sqrt(-D);
      cycubed = Complex(-q/2.0,sqrtD);

      A3 = hypot(cycubed.r,cycubed.i);  A = exp(log(max(MINREAL,A3))/3.0);
      theta3 = atan2(cycubed.i,cycubed.r); theta = theta3/3.0;
      cy1 = Complex(A*cos(theta),A*sin(theta));
      cy2 = Complex(A*cos(theta),-A*sin(theta));

      root[1] = Complex(cy1.r+cy2.r-scf[2]/3.0, cy1.i+cy2.i);
      root[2] = Cadd(Cmul(cy1,c3rtone),Cmul(cy2,Conjg(c3rtone)));
      root[2].r += - scf[2]/3.0;
      root[3] = Cadd(Cmul(cy1,Conjg(c3rtone)),Cmul(cy2,c3rtone));
      root[3].r += - scf[2]/3.0;
   }
}

#undef MINREAL
