/* ========================================================================
   This file contains various programs to interpolate functions.
 
      interpolate       - take derivative at interpolated point
      myspline/mysplint - use cubic spline
      getextremum       - find location and value of min/max

   ========================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"
#include "differentiate.h"
#include "alloc_space.h"
#include "myerror.h"
#include "constants.h"
#include "macros.h"

/* =========================================================
   Find interpolated point using linear fit
   =========================================================
 */
double interpolate0(xval, x1,x2,y1,y2)
   double xval,x1,x2,y1,y2;
{
  double d21;

  if ((d21=(x2-x1))==0.0) myerror("INTERPOLATE0: singularity encountered");

  return( (y2*(xval-x1) + y1*(x2-xval))/d21 );
}

/* =========================================================
   Find interpolated point  using cubic fit
   =========================================================
 */
double interpolate(xval, x1,x2,x3, y1,y2,y3)
   double xval,x1,x2,x3,y1,y2,y3;
{
  double d21,d32,d13,p11,p21,p22,p32,p33,p13,ipd123;
  double A,B,C;

  d21 = x2-x1;  d32 = x3-x2;  d13 = x1-x3;
  p11 = x1*x1;  p13 = x1*x3;
  p22 = x2*x2;  p21 = x2*x1;
  p33 = x3*x3;  p32 = x3*x2;
  ipd123 = 1.0/(d21*d32*d13);

  A = -ipd123*(y1*d32+y2*d13+y3*d21);
  B = ipd123*(y1*(p33-p22)+y2*(p11-p33)+y3*(p22-p11));
  C = -ipd123*(y1*p32*d32+y2*p13*d13+y3*p21*d21);

  return( (A*xval + B)*xval + C );
}

/* =====================================================================
   Given list (x[i],y[i]) and xval find y(xval) using linear 
   interpolation.  

   WARNING: Use of this procedure should ensure that y is monotonic in x
   =====================================================================
*/
void interpolate_list(x,y,n, xval,yval)
   int n;
   double xval;
   double *x,*y,*yval;
{
   int i;

   if (n<2) myerror("INTERPOLATE_LIST: n is too small");

   if (x[n]>x[1]) 
      if (xval<=x[1]) 
         *yval = interpolate0(xval,x[1],x[2],y[1],y[2]);
      else if (xval>=x[n]) 
         *yval = interpolate0(xval,x[n-1],x[n],y[n-1],y[n]);
      else {
         i=1;
         while (x[++i]<xval);
         *yval = interpolate0(xval,x[i-1],x[i],y[i-1],y[i]);
      }
   else  /* x[n]<=x[1]  (Note: x[n]=x[1] will give error) */
      if (xval>=x[1]) 
         *yval = interpolate0(xval,x[1],x[2],y[1],y[2]);
      else if (xval<=x[n]) 
         *yval = interpolate0(xval,x[n-1],x[n],y[n-1],y[n]);
      else {
         i=1;
         while (x[++i]>xval);
         *yval = interpolate0(xval,x[i-1],x[i],y[i-1],y[i]);
      }
}

/* ===================================================================
   FINDMAX:
  
   Find min/max point given three points.
   =================================================================== 
*/
int findmax(
   double x1, double x2, double x3, double y1, double y2, double y3,
   double *xm, double *ym

)
{
   double m21,m32,A,B,C;


   m21 = (y2-y1)/(x2-x1);
   m32 = (y3-y2)/(x3-x2);
   A = (m21-m32)/(x1-x3);
   B = m21 - (x1+x2)*A;
   C = y1 - (A*x1 + B)*x1;

   if (A==0) { /* points are linear */
     *xm=*ym=0.0;
     return(0);
   } else {
     *xm = -0.5*B/A;
     *ym = C-0.25*B*B/A;
     if (A>0) { /* minimum */ 
       return(-1);
     } else {   /* maximum */
       return(1);
     }
   }
}

/* =================================================================
   Find y2 containing data used for cubic spline 
   =================================================================
*/
void myspline(x,y,n, yp1,ypn,y2)
   int n;
   double yp1,ypn;
   double *x,*y,*y2;
{
	int i,k;
	double p,qn,sig,un,*u;

	u=dvector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_dvector(u,1,n-1);
}

/* ===================================================================
   Given data computed using myspline, calculate value of interpolated
   point using cubic spline.
   =================================================================== 
*/

void mysplint(xa,ya,y2a,n, x,y)
   int n;
   double x;
   double *xa,*ya,*y2a,*y;
{
	int klo,khi,k;
	double h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) myerror("Bad XA input to routine SPLINT");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


/* =================================================================
   INTERPOLATEFN: Nov 4, 2001
   =================================================================
*/
void interpolatefn(double *x,double *f,int n, 
                   double *x2,double *f2,int n2)
{
   int i;
   double xp1,xpN;
   double *fs;

   fs = dvector(1,n);

   xp1=intplderiv(x[1], x[1],x[2],x[3], f[1],f[2],f[3]);
   xpN=intplderiv(x[n], x[n],x[n-1],x[n-2], f[n],f[n-1],f[n-2]);

   myspline(x,f,n, xp1,xpN, fs);

   for (i=1;i<=n2;i++) mysplint(x,f,fs,n,x2[i],&f2[i]);
  
   free_dvector(fs,1,n);
}
