/* =================================================================

   This file contains routines to integrate functions of one or
      more variables.
 
      integrate - solve the ode between boundaries.  Returns number.
      intfdx    - simple trapezoidal sum method evaluated between
                  specified boundaries.  Returns number.
      intf      - simple trapezoidal sum method returning function
      integrateO1 - integrate 0..N grid using trapezoidal  method
   =================================================================
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"
#include "alloc_space.h"
#include "interpolate.h"
#include "myerror.h"
#include "constants.h"
#include "macros.h"

/* ================================================================
   INTEGRATE: (June 15, 1992)

   Program to integrate a function specified at evenly spaced points
   ================================================================
*/

#define TINY 1.0e-30
#define PGROW -0.20
#define PSHRNK -0.25
#define FCOR 0.06666666      /* 1.0/15.0 */
#define SAFETY 0.9
#define ERRCON 6.0e-4

#define Hstart 1.0e-3
#define Hmin 1.0e-20
#define MAXSTP 10000

void *int_derivs(xind, y,yp,npts, x,fx,fxspline)
   int npts;
   double xind,y;
   double *yp,*x,*fx,*fxspline;
{
   mysplint(x,fx,fxspline,npts,xind,yp);
}

void myrk4(y,dydx,x,h,yout, int_derivs, npts,z,fnz,fnzspline)
   int npts; 
   double y,dydx,x,h; 
   double *yout,*z,*fnz,*fnzspline;
   void (*int_derivs)();
{
   double xh,hh,h6,dym,dyt,yt;

   hh=h*0.5;
   h6=h/6.0;
   xh=x+hh;
   yt=y+hh*dydx;
   (*int_derivs)(xh,yt,&dyt,npts,z,fnz,fnzspline);
   yt=y+hh*dyt;
   (*int_derivs)(xh,yt,&dym,npts,z,fnz,fnzspline);

   yt=y+h*dym;
   dym += dyt;
   
   (*int_derivs)(x+h,yt,&dyt,npts,z,fnz,fnzspline);
   *yout=y+h6*(dydx+dyt+2.0*dym);
}

void myrkqc(y,dydx,x,htry,eps, yscal,hdid,hnext, int_derivs, 
      npts,z,fnz,fnzspline)
   int npts; 
   double htry,eps; 
   double *y,*dydx,*x,*yscal,*hdid,*hnext,*z,*fnz,*fnzspline;
   void (*int_derivs)();
{
   double xsav,hh,h,temp,errmax;
   double dysav,ysav,ytemp;

   xsav=(*x);
   
   ysav= *y;
   dysav= *dydx;
   
   h=htry;
   for (;;) {
      hh=0.5*h;
      myrk4(ysav,dysav,xsav,hh,&ytemp,int_derivs,npts,z,fnz,fnzspline);
      *x=xsav+hh;
       (*int_derivs)(*x,ytemp,dydx,npts,z,fnz,fnzspline);
      myrk4(ytemp,*dydx,*x,hh,y,int_derivs,npts,z,fnz,fnzspline);
      *x=xsav+h;
      if (*x == xsav) myerror("Step size too small in routine RKQC");
      myrk4(ysav,dysav,xsav,h,&ytemp,int_derivs,npts,z,fnz,fnzspline);
      errmax=0.0;
      
      ytemp= *y-ytemp;
      temp=dabs(ytemp/(*yscal));
      if (errmax < temp) errmax=temp;
      
      errmax /= eps;
      if (errmax <= 1.0) {
         *hdid=h;
         *hnext=(errmax > ERRCON ?
            SAFETY*h*exp(PGROW*log(errmax)) : 4.0*h);
         break;
      }
      h=SAFETY*h*exp(PSHRNK*log(errmax));
   }
   *y += ytemp*FCOR;
}

int kmax=0,kount=0;  /* defining declaration */
double *xp=0,*yp=0,dxsav=0;  /* defining declaration */
#define Epsintegrate 1.0e-6

double integrate(fx,x, npts, x1,x2)
   int npts;
   double x1,x2;
   double *fx,*x;
{
   int nok, nbad;
   int nstp;
   double integral;
   double xsav,xval,hnext,hdid,h;
   double yscal,y,dydx;
   double fxp1,fxpN;
   double *fxspline;

   fxspline =dvector(1,npts);
   fxp1 = (fx[2]-fx[1])/(x[2]-x[1]);
   fxpN = (fx[npts]-fx[npts-1])/(x[npts]-x[npts-1]);
    myspline(x,fx,npts,fxp1,fxpN,fxspline);

   xval=x1;
   h=(x2 > x1) ? Hstart : -Hstart;
   nok = nbad = kount = 0;
   y=0.0;
   if (kmax > 0) xsav=xval-dxsav*2.0;

   for (nstp=1;nstp<=MAXSTP;nstp++) {
      (*int_derivs)(xval,y,&dydx,npts,x,fx,fxspline);
      yscal=dabs(y)+dabs(dydx*h)+TINY;
      if (kmax > 0) {
         if (dabs(xval-xsav) > dabs(dxsav)) {
            if (kount < kmax-1) {
               xp[++kount]=xval;
               yp[kount]=y;
               xsav=xval;
            }
         }
      }
      if ((xval+h-x2)*(xval+h-x1) > 0.0) h=x2-xval;
      myrkqc(&y,&dydx,&xval,h,Epsintegrate,&yscal,&hdid,&hnext,int_derivs,
               npts,x,fx,fxspline); 
        if (hdid == h) ++nok; else ++nbad;
      if ((xval-x2)*(x2-x1) >= 0.0) {
         integral=y;
         if (kmax) {
            xp[++kount]=xval;
            yp[kount]=y;
         }
         free_dvector(fxspline,1,npts);
         return(integral);
      }
      if (dabs(hnext) <= Hmin) myerror("INTEGRATE: Step size too small");
      h=hnext;
   }
   myerror("INTEGRATE: Too many steps");
}

#undef MAXSTP
#undef Hmin
#undef Hstart
#undef TINY
#undef PGROW
#undef PSHRNK
#undef FCOR
#undef SAFETY
#undef ERRCON

/* ======================================================================
   rknstp,rknde, and myerny are also integrators which are suitable for
   differential equations of the form: y'' + f(x,y(x)) = 0.
   ======================================================================
*/

#define ERRCON 0.018      /* = (2/SAFETY)^(1/PGROW) */
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define TINY   1.e-20
#define MACHTINY 0.000000000000000222

/* variables for myerny routine */
#define Hstart 0.001
#define Hmin 1.0e-20
#define Hmax 10.0
#define MAXINTSTP 50000


void rknstp (y,yp,d2ydx2, nvar,x, h,eps,hdid,hnext, derivs, 
      odeparams,flags,zind,U,U2,Upp,Upp2,nsqr,nsqr2,nind)
   int nvar, nind;
   int *flags;
   double h,eps;
   double *y,*yp,*d2ydx2,*x,*hdid,*hnext,*odeparams;
   double zind[],U[],U2[],Upp[],Upp2[],nsqr[],nsqr2[];
   void (*derivs)();
/*
-------------------------- prologue -----------------------------

    This subroutine is the step integrator and can either be used
  by itself or with the inner driver rknde. the former way permits the
  greater control over the integration but requires more effort on
  behalf of the user.
     rknstp attempts a step of the local problem given y and y'' at x.
  if the step is accepted the final values will be returned in y and yp.
  
  Two machine dependent numbers are used in rknstp. the
  unit roundoff, u, and a very small positive machine number, TINY. 
  The values for both are obtained by calling the subroutine machne.
*/
{
   int i;
   double c2,c3,c4,c5,a21,a31,a32,a41,a42,a43,a51,a52;
   double a53,a54,a61,a63,a64,a65,bp1,bp3,bp4,bp5,bp6,d1,d3,d4,d5,dp1;
   double dp3,dp4,dp5,dp6;
   double c2h,c3h,c4h,c5h;
   double erry,erryp,errmax,hsqd;
   double tmpx;
   double *tmpy, *tmpyp, **tmpypp;

   tmpy = dvector(1,nvar);
   tmpyp = dvector(1,nvar);
   tmpypp = dmatrix2(1,6,1,nvar);

/*  runge-kutta nystrom coefficients */

   c2 = 2.0/15.0;     c3 = 4.0/15.0;
   c4 = 11.0/16.0;    c5 = 6.0/7.0;
   
   a21 = 2.0/225.0;            a31 = 8.0/675.0;
   a32 = 16.0/675.0;           a41 = 315931.0/2097152.0;
   a42 = -246235.0/1048576.0;   a43 = 672155.0/2097152.0;
   a51 = -1040407.0/3565485.0;  a52 = 8.0/9.0;
   a53 = -748807.0/2182509.0;   a54 = 40819712.0/360113985.0;
   a61 = 85.0/1056.0;          a63 = 29475.0/100192.0;
   a64 = 7168.0/63327.0;       a65 = 343.0/28272.0;

   bp1 = 85.0/1056.0;          bp3 = 442125.0/1102112.0;
   bp4 = 114688.0/316635.0;    bp5 = 2401.0/28272.0;
   bp6 = 47.0/660.0;
     
   d1 = -7.0/1056.0;           d3 = 17325.0/1102112.0;
   d4 = -6720.0/316635.0;      d5 = 343.0/28272.0;
     
   dp1 = 47.0/3168.0;          dp3 = -52875.0/1102112.0;
   dp4 = 48128.0/316635.0;     dp5 = -16121.0/84816.0;
   dp6 = 47.0/660.0;

   /*  MAIN LOOP - loop until successful step or error condition */

   for (;;) { /* perform derivative evaluations */
      if (h>Hmax) 
         myerror( "Step size is too large in RKNSTP" );

      hsqd  = h*h;

      /* tmpypp1 */

      (*derivs)(*x,y,tmpypp[1],odeparams,flags,zind,U,U2,Upp,Upp2,nsqr,nsqr2,nind);

      /* tmpypp2 */

      c2h  = c2*h;
      tmpx = *x + c2h;
      for (i=1;i<=nvar;i++) 
         tmpy[i] = y[i] + c2h*yp[i] + hsqd*a21*tmpypp[1][i];
      (*derivs)(tmpx,tmpy,tmpypp[2],odeparams,flags,zind,U,U2,Upp,Upp2,nsqr,nsqr2,nind);

      /* tmpypp3 */

      c3h  = c3*h;
      tmpx = *x + c3h;
      for (i=1;i<=nvar;i++) 
         tmpy[i] = y[i] + c3h*yp[i] + hsqd*(a31*tmpypp[i][1]+a32*tmpypp[i][2]);
      (*derivs)(tmpx,tmpy,tmpypp[3],odeparams,flags,zind,U,U2,Upp,Upp2,nsqr,nsqr2,nind);

      /* tmpypp4 */

      c4h  = c4*h;
      tmpx = *x + c4h;
      for (i=1;i<=nvar;i++) 
         tmpy[i] = y[i] + c4h*yp[i] + 
                   hsqd*(a41*tmpypp[1][i]+a42*tmpypp[2][i]+a43*tmpypp[3][i]);
     (*derivs)(tmpx,tmpy,tmpypp[4],odeparams,flags,zind,U,U2,Upp,Upp2,nsqr,nsqr2,nind);

      /* tmpypp5 */

      c5h  = c5*h;
      tmpx = *x + c5h;
      for (i=1;i<=nvar;i++) 
         tmpy[i] = y[i] + c5h*yp[i] + 
                   hsqd*( a54*tmpypp[4][i]+a51*tmpypp[1][i]+ 
                          a53*tmpypp[3][i]+a52*tmpypp[2][i] );
      (*derivs)(tmpx,tmpy,tmpypp[5],odeparams,flags,zind,U,U2,Upp,Upp2,nsqr,nsqr2,nind);

      /* tmpypp6 */

      tmpx = *x + h;
      for (i=1;i<=nvar;i++) 
         tmpy[i] = y[i] + h*yp[i] + 
                   hsqd*( a65*tmpypp[5][i]+a61*tmpypp[1][i]+ 
                          a64*tmpypp[4][i]+a63*tmpypp[3][i] );
      (*derivs)(tmpx,tmpy,tmpypp[6],odeparams,flags,zind,U,U2,Upp,Upp2,nsqr,nsqr2,nind);

      if (*x == tmpx)  /* check if step size is too small */ 
         myerror("Step size too small in RKNSTP");

      /* calculate the norm of the local error estimate */

      errmax = 0.0;
      for (i=1;i<=nvar;i++) { 
         erry  = hsqd*(dabs(d1*tmpypp[1][i]+d3*tmpypp[3][i]+
                           d4*tmpypp[4][i]+d5*tmpypp[5][i]));
         erryp = (dabs(h))*(dabs(dp1*tmpypp[1][i]+dp3*tmpypp[3][i]+
                              dp4*tmpypp[4][i]+dp5*tmpypp[5][i]+
                              dp6*tmpypp[6][i]));
         errmax = max( errmax, max( erry, erryp ));
      }

      /* if the step is to be accepted update x,y,y', select a new
         stepsize and exit.  The step is rejected if the size is
         too small.  Otherwise select a smaller stepsize and continue
      */

      errmax /= eps;
      if (errmax <= 1.0) {  /* step accepted */
         *x += h;
         *hdid = h;
         *hnext = (errmax>ERRCON ? (SAFETY*h*exp(PGROW*log(errmax))) : 2.0*h);
         break;
      } 
      h = SAFETY*h*exp(PSHRNK*log(errmax));
   }
   for (i=1;i<=nvar;i++) {
       y[i] = tmpy[i];
       yp[i] += h*( bp1*tmpypp[1][i] + bp3*tmpypp[3][i] + 
                    bp4*tmpypp[4][i] + bp5*tmpypp[5][i] + bp6*tmpypp[6][i]);
   }
   free_dvector(tmpy,1,nvar);
   free_dvector(tmpyp,1,nvar);
   free_dmatrix2(tmpypp,1,6,1,nvar);
}


void rknde (ystart,ypstart, nvar,xin,xout, eps,nok,nbad, derivs, 
      odeparams,flags, zind,U,U2,Upp,Upp2,nsqr,nsqr2, nind,
      printflg,normphi,efn,tnstp)
   int nvar,printflg;
   int *nok,*nbad,*flags,*tnstp;
   double xin,xout,eps;
   double *ystart,*ypstart,*odeparams,*normphi;
   double zind[],U[],U2[],Upp[],Upp2[],nsqr[],nsqr2[];
   double **efn;
   void (*derivs)();

{
/* ---------------------------- prologue -------------------------------
 
        this subroutine is the inner driver for the integrator erny. it
     is intended to be used with the outer driver which merely checks
     the input values and allocates storage. rknde handles the
     information for the step integrator rknstp and the interpolation
     routine intrp together with some detection of improper use of
     erny. a more detailed description of the use of rknde is given in
     the prologue for the outer driver (erny).
*/
   int i,nstp;
   double h,hdid,hnext,x;
   double *y,*yp,*d2ydx2;
     
   y = dvector(1,nvar);
   yp = dvector(1,nvar);
   d2ydx2 = dvector(1,nvar);

   h = (xout>xin) ? Hstart : - Hstart;
   x = xin;
   for (i=1;i<=nvar;i++) {
      y[i] = ystart[i];
      yp[i] = ypstart[i];
   }
   if (printflg!=0) {
      efn[0][0] = x;
      efn[0][1] = normphi[1]*y[1];
      efn[0][2] = normphi[2]*yp[1];
      efn[0][3] = normphi[3]*y[2];
      efn[0][4] = normphi[4]*yp[2];
   }

   for (nstp=1;nstp<=MAXINTSTP;nstp++) { /* loop over integration steps */
      (*derivs)(x,y,d2ydx2,odeparams,flags,zind,U,U2,Upp,Upp2,nsqr,nsqr2,nind);

      /* If the integration is past the output point, reduce step size.
         Only check for this after first pass since initially h could
         equal zero */

      if ( (x+h-xout)*(x+h-xin) >= 0.0 ) h = xout-x;

      /*  take a step */

      rknstp(y,yp,d2ydx2,nvar,&x, h,eps,&hdid,&hnext, derivs,odeparams,flags,
             zind,U,U2,Upp,Upp2,nsqr,nsqr2,nind);
 
      /* if required, print (normalized) results */
      if (printflg!=0) {
         efn[nstp][0] = x;
         efn[nstp][1] = normphi[1]*y[1];
         efn[nstp][2] = normphi[2]*yp[1];
         efn[nstp][3] = normphi[3]*y[2];
         efn[nstp][4] = normphi[4]*yp[2];
      }
      if (printflg%10 ==2) 
         printf("%8.6lf\t%8.6lf\t%8.6lf\t%8.6lf\t%8.6lf\n",
             efn[nstp][0],efn[nstp][1],efn[nstp][2],efn[nstp][3],efn[nstp][4]);

      /* check if step size used equals size set */

      if (hdid == h) ++(*nok); else ++(*nbad);

      /* check if finished */
      if ( (dabs(x-xout)) < eps) {
         for (i=1;i<=nvar;i++) {
            ystart[i] = y[i];
            ypstart[i] = yp[i];
         }
         free_dvector(d2ydx2,1,nvar);
         free_dvector(y,1,nvar);
         free_dvector(yp,1,nvar);
         *tnstp = nstp;
         return;
      }
      
      /* check if step size is too small */
      if (dabs(hnext) <= Hmin) 
         myerror("Step size too small in RKNDE");
      h = hnext;
   }
   myerror( "Integration required more steps than MAXINTSTP, in RKNDE\n" );
}

void myerny(phi,nvar,xin,xout, eps,nok,nbad, derivs, 
      odeparams,flags, zind,U,U2,Upp,Upp2,nsqr,nsqr2, nind,
      printflg,normphi,efn,tnstp)
   int nvar,nind,printflg;
   int *nok,*nbad,*flags,*tnstp;
   double xin,xout,eps;
   double *phi,*odeparams,*normphi;
   double zind[],U[],U2[],Upp[],Upp2[],nsqr[],nsqr2[];
   double **efn;
   void (*derivs)();
{
/*
 --------------------------- prologue --------------------------------

      The integrator erny finds approximate values for the solution and
   derivative of the second order initial value problem:

             y'' = f(x,y),  y(x0) = y0, y'(x0) = y'0
             where y,f are n-vectors, ' denotes d/dx

   at the output point x = xout. upon returning with these values the
   integration can be continued merely by changing xout and recalling
   erny. in this way y and y' can be approximated at a sequence of
   output points. the code attempts to keep the error in y and y'
   proportional to a prescribed epserance and is intended for smooth
   problems. if f is expensive to evaluate or high accuracy is required
   a variable-order variable-step adams or variable-step high order
   runge-kutta nystrom may be more appropriate.
      The integration is based on a variable-step fourth and fifth
   order explicit runge-kutta nystrom pair developed by Bettis and Horn.
   Unless specified differently by the user,
   the code will integrate past the output point and use interpolation
   to find the values of y and y' at xout. the user can control several
   aspects of the integration including the type of error control and
   the maximum number of derivative evaluations that are permitted. 
*/
   int i;
   double *ystart, *ypstart;

/* check the input values and exit if the values are invalid.
         ------ nvar > 0, eps > 0 --------
*/

   ystart = dvector(1,nvar);
   ypstart = dvector(1,nvar);
   for (i=1;i<=nvar;i++) {
      ystart[i] = phi[2*i-1];
      ypstart[i] = phi[2*i];
   }

   /*  ---- attempt the integration ---- */

   rknde(ystart, ypstart, nvar, xin, xout, eps, nok, nbad,
           derivs, odeparams, flags, zind, U,U2,Upp,Upp2,nsqr,nsqr2,nind,
           printflg, normphi, efn,tnstp);

   /* clean up variables */
   for (i=1;i<=nvar;i++) {
      phi[2*i-1] = ystart[i];
      phi[2*i] = ypstart[i];
   }
   free_dvector(ystart,1,nvar);
   free_dvector(ypstart,1,nvar);
}

#undef ERRCON
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef TINY
#undef MACHTINY
#undef Hstart
#undef Hmin 
#undef Hmax
#undef MAXINTSTP

/* ==================================================================
   END OF MYERNY AND ASSOCIATED SUBROUTINES
   ==================================================================
*/


/* ==================================================================
   INTFDX: (June 17, 1992)

   Program to integrate a function specified at evenly spaced points
   Integration proceeds by simple trapezoidal method.
   ==================================================================
*/

double intfdx(x,f, npts, x1,x2)
   int npts;
   double x1,x2;
   double *f,*x;
{
   int i,indl,indh;
   double sum;

   indl=1;
   while (x[indl]<x1) indl++;
   indh=npts;
   while (x[indh]>x2) indh--;

   sum=0.0;
   for (i=indl;i<indh;i++) 
      sum += 0.5*(f[i+1]+f[i]) * (x[i+1]-x[i]);

   return(sum);
}

/* =============================================================
   INTF:

   Simple routine to integrate f[i](x[i]) using trapezoidal sums
   =============================================================
 */
void intf(fn,x, nind,ifn)
   int nind;
   double *fn,*x,*ifn;
{
  int i;

  ifn[1]=0.0;
  for (i=2;i<=nind;i++) 
     ifn[i] = ifn[i-1]+ 0.5*(fn[i]+fn[i-1])*(x[i]-x[i-1]);
}

/* =========================================================================
   INTEGRATEO0: Oct 10/93

   Find integral of a real function given initial value and length of
   domain.  f is  defined for indices from 0 to N at equally spaced
   intervals of size dx.  Integrate forward assigning if_i = f_i.
   =========================================================================
*/
void integrateO0(f, ox,nx, dx,if0, irf)
   int ox,nx;
   double dx,if0;
   double *f,*irf;
{
   int i;

   irf[ox] = if0+dx*f[ox];
   for (i=ox+1;i<=nx;i++) irf[i] = irf[i-1] + dx*f[i];
}

/* =========================================================================
   INTEGRATEO1: Sept 28/93

   Find integral of a real function given initial value and length of
   domain.  f is  defined for indices from 0 to N at equally spaced
   intervals of size dx.

   If isign==1 integrate forward assigning if_i = 0.5(f_i+1 + f_i)
   If isign==-1 integrate forward assigning if_i = 0.5 (f_i + f_i-1)
   The initial value of the integral is taken to be if0.
   =========================================================================
*/
void integrateO1(f, ox,nx, dx,isign,if0, irf)
   int ox,nx,isign;
   double dx,if0;
   double *f,*irf;
{
   int i;
   double hfdx;

   if (iabs(isign)!=1)
      myerror("INTEGRATEO1: isign must equal +/-1");

   hfdx = 0.5*dx;
   if (isign==1) {
      irf[ox] = if0+hfdx*(f[ox]+f[ox+1]);
      for (i=ox+1;i<nx;i++) 
         irf[i] = irf[i-1] + hfdx*(f[i]+f[i+1]);
      irf[nx] = 2*irf[nx-1] - irf[nx-2];
      if (ox==1)
         irf[0] = 2*irf[1] - irf[2];
   } else if (isign==-1) {
      irf[ox]=if0;
      for (i=ox+1;i<=nx;i++) 
         irf[i] = irf[i-1] + hfdx*(f[i]+f[i-1]);
   }
}

/* =========================================================================
   INTEGRATEO1_gen: Sept 28/93

   Find integral of a real function given initial value and length of
   domain.  x,f is  defined for indices from 1 to N.

   if isign==1 integrate forward  if[0] = if0 
   if isign==-1 integrate backward  if[n] = if0
   =========================================================================
*/
void integrateO1_gen(f,x,ox,nx, isign,if0, irf)
   int ox,nx,isign;
   double if0;
   double *f,*x,*irf;
{
   int i;
   double hfdx;

   if (iabs(isign)!=1)
      myerror("INTEGRATEO1_GEN: isign must equal +/- 1");

   if (isign==1) {
      irf[ox]=if0;
      for (i=ox+1;i<=nx;i++) {
         hfdx = 0.5*(x[i]-x[i-1]);
         irf[i] = irf[i-1] + hfdx*(f[i]+f[i-1]);
      }
   } else if (isign==-1) {
      irf[nx]=if0;
      for (i=nx-1;i>=ox;i--) {
         hfdx = 0.5*(x[i+1]-x[i]);
         irf[i] = irf[i+1] + hfdx*(f[i]+f[i+1]);
      }
   }
}

/* =========================================================================
   CINTEGRATEO0: Oct 11/93

   Find integral of a real function given initial value and length of
   domain.  f is  defined for indices from 0 to N at equally spaced
   intervals of size dx.  Integrate forward assigning if_i = f_i.
   =========================================================================
*/
void CintegrateO0(f,ox,nx, dx,isign,if0, icf)
   int ox,nx,isign;
   double dx;
   dcomplex if0;
   dcomplex *f,*icf;
{
   int i;

   if (iabs(isign) != 1)
      myerror("CINTEGRATEO1: isign must equal +/- 1");

   if (isign==1) {
      icf[ox].r = if0.r;
      icf[ox].i = if0.i;
      for (i=ox+1;i<=nx;i++) {
         icf[i].r = icf[i-1].r + dx*f[i-1].r;
         icf[i].i = icf[i-1].i + dx*f[i-1].i;
      }
   } else {
      icf[ox].r = if0.r+dx*f[ox].r;
      icf[ox].i = if0.i+dx*f[ox].i;
      for (i=ox+1;i<=nx;i++) {
         icf[i].r = icf[i-1].r + dx*f[i].r;
         icf[i].i = icf[i-1].i + dx*f[i].i;
      }
   }
}

/* =========================================================================
   CINTEGRATEO0_gen: Nov 20/93

   Find integral of a complex function given initial value and length of
   domain.  x,cf is  defined for indices from 1 to N.

   if isign==1 integrate forward  icf[0] = icf0 
   if isign==-1 integrate backward  icf[n] = icf0
   =========================================================================
*/
void CintegrateO0_gen(cf,x,ox,nx, isign,icf0, icf)
   int ox,nx,isign;
   double *x;
   dcomplex icf0;
   dcomplex *cf,*icf;
{
   int i;
   double dx;

   if (iabs(isign) != 1)
      myerror("CINTEGRATEO1_GEN: isign must equal +/- 1");

   if (isign==1) {
      icf[ox]=icf0;
      for (i=ox+1;i<=nx;i++) {
         dx = x[i]-x[i-1];
         icf[i].r = icf[i-1].r + dx*cf[i-1].r;
         icf[i].i = icf[i-1].i + dx*cf[i-1].i;
      }
   } else if (isign==-1) {
      icf[nx]=icf0;
      for (i=nx-1;i>=ox;i--) {
         dx = x[i-1]-x[i];
         icf[i].r = icf[i-1].r + dx*cf[i].r;
         icf[i].i = icf[i-1].i + dx*cf[i].i;
      }
   }
}

/* =========================================================================
   CINTEGRATEO1: SEPT 28/93

   Find integral of a complex function given initial value and length of
   domain.  cf is  defined for indices from 1 to N at equally spaced
   intervals of length dx.

   If isign==1 integrate forward assigning icf_i = 0.5*(cf_i + cf_i+1)
   If isign==-1 integrate forward assigning icf_i = 0.5*(cf_i + cf_i-1)
   Initial value is taken to be icf0.
   =========================================================================
*/
void CintegrateO1(cf, ox,nx, dx,isign,icf0, icf)
   int ox,nx,isign;
   double dx;
   dcomplex icf0;
   dcomplex *cf,*icf;
{
   int i;
   double hfdx;

   if (iabs(isign) != 1)
      myerror("CINTEGRATEO1: isign must equal +/- 1");

   hfdx = 0.5*dx;
   if (isign==1) {
      icf[ox].r=icf0.r + hfdx*(cf[ox].r+cf[ox+1].r);
      icf[ox].i=icf0.i + hfdx*(cf[ox].i+cf[ox+1].i);
      for (i=ox+1;i<nx;i++) {
         icf[i].r = icf[i-1].r + hfdx*(cf[i].r+cf[i+1].r);
         icf[i].i = icf[i-1].i + hfdx*(cf[i].i+cf[i+1].i);
      }
      icf[nx].r = 2*icf[nx-1].r - icf[nx-2].r;
      icf[nx].i = 2*icf[nx-1].i - icf[nx-2].i;
      if (ox==1) {
         icf[0].r = 2*icf[1].r - icf[2].r;
         icf[0].i = 2*icf[1].i - icf[2].i;
      }
   } else if (isign==-1) {
      icf[ox]=icf0;
      for (i=ox+1;i<=nx;i++) {
         icf[i].r = icf[i-1].r + hfdx*(cf[i].r+cf[i-1].r);
         icf[i].i = icf[i-1].i + hfdx*(cf[i].i+cf[i-1].i);
      }
   }
}

/* =========================================================================
   CINTEGRATEO1_gen: SEPT 28/93

   Find integral of a complex function given initial value and length of
   domain.  x,cf is  defined for indices from 1 to N.

   if isign==1 integrate forward  icf[0] = icf0 
   if isign==-1 integrate backward  icf[n] = icf0
   =========================================================================
*/
void CintegrateO1_gen(cf,x, ox,nx, isign,icf0, icf)
   int ox,nx,isign;
   double *x;
   dcomplex icf0;
   dcomplex *cf,*icf;
{
   int i;
   double hfdx;

   if (iabs(isign) != 1)
      myerror("CINTEGRATEO1_GEN: isign must equal +/- 1");

   if (isign==1) {
      icf[ox]=icf0;
      for (i=ox+1;i<=nx;i++) {
         hfdx = 0.5*(x[i]-x[i-1]);
         icf[i].r = icf[i-1].r + hfdx*(cf[i].r+cf[i-1].r);
         icf[i].i = icf[i-1].i + hfdx*(cf[i].i+cf[i-1].i);
      }
   } else if (isign==-1) {
      icf[nx]=icf0;
      for (i=nx-1;i>=ox;i--) {
         hfdx = 0.5*(x[i+1]-x[i]);
         icf[i].r = icf[i+1].r + hfdx*(cf[i].r+cf[i+1].r);
         icf[i].i = icf[i+1].i + hfdx*(cf[i].i+cf[i+1].i);
      }
   }
}



/* =========================================================================
   ZINTEGRATEO0: Oct 11/93

   Find integral of a real function given initial value and length of
   domain.  f is  defined for indices from 0 to N at equally spaced
   intervals of size dx.  Integrate forward assigning if_i = f_i.
   =========================================================================
*/
void ZintegrateO0(f,ox,nx, dx,isign,if0, icf)
   int ox,nx,isign;
   double dx;
   zomplex if0;
   zomplex *f,*icf;
{
   int i;

   if (iabs(isign) != 1)
      myerror("ZINTEGRATEO1: isign must equal +/- 1");

   if (isign==1) {
      icf[ox].re = if0.re;
      icf[ox].im = if0.im;
      for (i=ox+1;i<=nx;i++) {
         icf[i].re = icf[i-1].re + dx*f[i-1].re;
         icf[i].im = icf[i-1].im + dx*f[i-1].im;
      }
   } else {
      icf[ox].re = if0.re+dx*f[ox].re;
      icf[ox].im = if0.im+dx*f[ox].im;
      for (i=ox+1;i<=nx;i++) {
         icf[i].re = icf[i-1].re + dx*f[i].re;
         icf[i].im = icf[i-1].im + dx*f[i].im;
      }
   }
}

/* =========================================================================
   ZINTEGRATEO0_gen: Nov 20/93

   Find integral of a complex function given initial value and length of
   domain.  x,cf is  defined for indices from 1 to N.

   if isign==1 integrate forward  icf[0] = icf0 
   if isign==-1 integrate backward  icf[n] = icf0
   =========================================================================
*/
void ZintegrateO0_gen(cf,x,ox,nx, isign,icf0, icf)
   int ox,nx,isign;
   double *x;
   zomplex icf0;
   zomplex *cf,*icf;
{
   int i;
   double dx;

   if (iabs(isign) != 1)
      myerror("ZINTEGRATEO1_GEN: isign must equal +/- 1");

   if (isign==1) {
      icf[ox]=icf0;
      for (i=ox+1;i<=nx;i++) {
         dx = x[i]-x[i-1];
         icf[i].re = icf[i-1].re + dx*cf[i-1].re;
         icf[i].im = icf[i-1].im + dx*cf[i-1].im;
      }
   } else if (isign==-1) {
      icf[nx]=icf0;
      for (i=nx-1;i>=ox;i--) {
         dx = x[i-1]-x[i];
         icf[i].re = icf[i-1].re + dx*cf[i].re;
         icf[i].im = icf[i-1].im + dx*cf[i].im;
      }
   }
}

/* =========================================================================
   ZINTEGRATEO1: SEPT 28/93

   Find integral of a complex function given initial value and length of
   domain.  cf is  defined for indices from 1 to N at equally spaced
   intervals of length dx.

   If isign==1 integrate forward assigning icf_i = 0.5*(cf_i + cf_i+1)
   If isign==-1 integrate forward assigning icf_i = 0.5*(cf_i + cf_i-1)
   Initial value is taken to be icf0.
   =========================================================================
*/
void ZintegrateO1(cf, ox,nx, dx,isign,icf0, icf)
   int ox,nx,isign;
   double dx;
   zomplex icf0;
   zomplex *cf,*icf;
{
   int i;
   double hfdx;

   if (iabs(isign) != 1)
      myerror("ZINTEGRATEO1: isign must equal +/- 1");

   hfdx = 0.5*dx;
   if (isign==1) {
      icf[ox].re=icf0.re + hfdx*(cf[ox].re+cf[ox+1].re);
      icf[ox].im=icf0.im + hfdx*(cf[ox].im+cf[ox+1].im);
      for (i=ox+1;i<nx;i++) {
         icf[i].re = icf[i-1].re + hfdx*(cf[i].re+cf[i+1].re);
         icf[i].im = icf[i-1].im + hfdx*(cf[i].im+cf[i+1].im);
      }
      icf[nx].re = 2*icf[nx-1].re - icf[nx-2].re;
      icf[nx].im = 2*icf[nx-1].im - icf[nx-2].im;
      if (ox==1) {
         icf[0].re = 2*icf[1].re - icf[2].re;
         icf[0].im = 2*icf[1].im - icf[2].im;
      }
   } else if (isign==-1) {
      icf[ox]=icf0;
      for (i=ox+1;i<=nx;i++) {
         icf[i].re = icf[i-1].re + hfdx*(cf[i].re+cf[i-1].re);
         icf[i].im = icf[i-1].im + hfdx*(cf[i].im+cf[i-1].im);
      }
   }
}

/* =========================================================================
   ZINTEGRATEO1_gen: SEPT 28/93

   Find integral of a complex function given initial value and length of
   domain.  x,cf is  defined for indices from 1 to N.

   if isign==1 integrate forward  icf[0] = icf0 
   if isign==-1 integrate backward  icf[n] = icf0
   =========================================================================
*/
void ZintegrateO1_gen(cf,x, ox,nx, isign,icf0, icf)
   int ox,nx,isign;
   double *x;
   zomplex icf0;
   zomplex *cf,*icf;
{
   int i;
   double hfdx;

   if (iabs(isign) != 1)
      myerror("ZINTEGRATEO1_GEN: isign must equal +/- 1");

   if (isign==1) {
      icf[ox]=icf0;
      for (i=ox+1;i<=nx;i++) {
         hfdx = 0.5*(x[i]-x[i-1]);
         icf[i].re = icf[i-1].re + hfdx*(cf[i].re+cf[i-1].re);
         icf[i].im = icf[i-1].im + hfdx*(cf[i].im+cf[i-1].im);
      }
   } else if (isign==-1) {
      icf[nx]=icf0;
      for (i=nx-1;i>=ox;i--) {
         hfdx = 0.5*(x[i+1]-x[i]);
         icf[i].re = icf[i+1].re + hfdx*(cf[i].re+cf[i+1].re);
         icf[i].im = icf[i+1].im + hfdx*(cf[i].im+cf[i+1].im);
      }
   }
}
