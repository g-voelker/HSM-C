#include <stdio.h>
#include <math.h>
#include "macros.h"
#include "constants.h"

/* ===========================================================================
   HeavysideAverage (Jun 17/96):

   Routine to take function with finely (irregularly) spaced abscissae, and
   return function defined with coarse, regularly spaced abscissae.
   This function works simply by averaging the original function over
   each interval defined by the coarse grid.
   ===========================================================================
*/
void HeavysideAverage(
   double *xraw, double *fxraw, int nraw,
   double lx,double range,
   double *x, double *fx, int nn)
{
  int i,nsum,n;
  double sum,hfrng;

  hfrng = 0.5*range;

  for (n=1;n<=nn;n++) {
     nsum=0;
     sum=0.0;
     for (i=1;i<=nraw;i++) {
        if (xraw[i]>=x[n]-hfrng && xraw[i]<=x[n]+hfrng) {
           sum += fxraw[i];
           nsum++;
        }
     }
     if (nsum>0) {
        fx[n] = sum/((double) nsum);
     } else {
        fx[n] = (n==1 ? 0.0 : fx[n-1]);
     }
  }
}

/* ===========================================================================
   GaussianAverage (Jun 17/96):
                   (Modified Apr 25/00)

   Routine to take function with finely (irregularly) spaced abscissae, and
   return function defined with coarse, regularly spaced abscissae.
   This function works simply by averaging the original function over
   each interval defined by the coarse grid.
   ===========================================================================
*/
void GaussianAverage(
   double *xraw, double *fxraw, int nraw,
   double dx,double range,
   double *x, double *fx, double *sd, int nn)
{
  int i,n;
  double sum,wt,wtsum,rngsqr,sumsum;

  rngsqr = range*range;

  for (n=1;n<=nn;n++) {
     sumsum=sum=wtsum=0.0;
     for (i=1;i<=nraw;i++) {
        wt = exp(-0.5*(x[n]-xraw[i])*(x[n]-xraw[i])/rngsqr);
        sum += wt*fxraw[i];
        sumsum += wt*fxraw[i]*fxraw[i];
        wtsum += wt;
     }
     fx[n] = sum/wtsum;
     sd[n] = sqrt(sumsum/wtsum - fx[n]*fx[n]);
  }
}

/* ===========================================================================
   BestFitLineSmooth (Jun 17/96):

   Routine to take function with finely (irregularly) spaced abscissae, and
   return function defined with coarse, regularly spaced abscissae.
   This function works by finding the best fit line over each interval
   defined by the coarse grid and evaluating the line at the midpoint.
   ===========================================================================
*/
void BestFitLineSmooth(
   double *xraw, double *fxraw, int nraw,
   double lx,double range,
   double *x, double *fx, int nn)
{
  int i,nsum,n;
  double xsum,zsum,xxsum,xzsum,m,b,hfrng;

  hfrng = 0.5*range;

  for (n=1;n<=nn;n++) {
     nsum = 0;
     zsum = xsum = xzsum = xxsum = 0.0;
     for (i=1;i<=nraw;i++) {
        if (xraw[i]>=x[n]-hfrng && xraw[i]<=x[n]+hfrng) {
           xsum += xraw[i];
           zsum += fxraw[i];
           xzsum += xraw[i]*fxraw[i];
           xxsum += xraw[i]*xraw[i];
           nsum++;
        }
     }

     if (nsum>0) {
        m = (nsum*xzsum - zsum*xsum)/(nsum*xxsum - xsum*xsum);
        b = (zsum - m*xsum)/((double) nsum);
        fx[n] = m*x[n] + b; 
     } else {
        fx[n] = (n==1 ? 0.0 : fx[n-1]);
     }
  }
}

#define MINREAL 1.0e-8

/* ===========================================================================
   BestFitLineFilter (Jun 17/96):

   Routine to take function with finely (irregularly) spaced abscissae, and
   return function defined with coarse, regularly spaced abscissae.
   This function works by finding the best fit line over each interval
   defined by the coarse grid, then throwing away points lying well outside
   the calculated std. dev., then re-evaluating the best fit line of
   this filtered data.  The value of this line at midpoint in interval is
   returned for each grid point.
   ===========================================================================
*/
void BestFitLineFilter(
   double *xraw, double *fxraw, int nraw,
   double lx,double range, double nsdm, double nsdb,
   double *x, double *fx, int nn)
{
  int filtered;
  int i,nsum,n,nfiltered;
  double denom,hfrng;
  double xsum,zsum,xxsum,xzsum,m,b;
  double meanerrsum,meanerr,merr,berr,minfx,maxfx;

  hfrng = 0.5*range;

  for (n=1;n<=nn;n++) {
     filtered = FALSE;
     nsum = 0;
     zsum = xsum = xzsum = xxsum = 0.0;
     for (i=1;i<=nraw;i++) {
        if (xraw[i]>=x[n]-hfrng && xraw[i]<=x[n]+hfrng) {
           xsum += xraw[i];
           zsum += fxraw[i];
           xzsum += xraw[i]*fxraw[i];
           xxsum += xraw[i]*xraw[i];
           nsum++;
        }
     }

     if (nsum>0) { /* find slope and intercept of best fit line */
        denom = nsum*xxsum - xsum*xsum;
        if (dabs(denom)<MINREAL) {  /* all values at one x point */
           fx[n] = zsum/((double) nsum);
        } else {
           m = (nsum*xzsum - zsum*xsum)/denom;
           b = (zsum - m*xsum)/((double) nsum);

           if (nsum<=2) {  /* Return line for unfiltered data */
              fx[n] = m*x[n] + b;
           } else {  /* Calculate errors for best fit line */
              filtered = TRUE;
              meanerrsum = 0;
              for (i=1;i<=nraw;i++)
                 if (xraw[i]>=x[n]-hfrng && xraw[i]<x[n]+hfrng) 
                   meanerrsum += (fxraw[i]-b-m*xraw[i])*(fxraw[i]-b-m*xraw[i]);

              meanerr = sqrt( meanerrsum/((double) (nsum-2)) );
              merr = meanerr*sqrt( xxsum/(nsum*xxsum - xsum*xsum) );
              berr = meanerr*sqrt( nsum/(nsum*xxsum - xsum*xsum) );
           }
        }
     }

     if (filtered) {
        /* now recalculate best fit line for points with std. dev. */
        nsum = 0;
        zsum = xsum = xzsum = xxsum = 0.0;
        for (i=1;i<=nraw;i++) {
           minfx = (m-nsdm*merr)*xraw[i] + b-nsdb*berr;
           maxfx = (m+nsdm*merr)*xraw[i] + b+nsdb*berr;
           if (xraw[i]>=x[n]-hfrng && xraw[i]<=x[n]+hfrng &&
               fxraw[i]>minfx && fxraw[i]<maxfx ) {
              xsum += xraw[i];
              zsum += fxraw[i];
              xzsum += xraw[i]*fxraw[i];
              xxsum += xraw[i]*xraw[i];
              nsum++;
           }
        }

        if (nsum==0) { /* return value for fit line to unfiltered data */
           fx[n] = m*x[n] + b;
        } else {       /* recalculate m and b for fit line to filtered data */
           denom = nsum*xxsum - xsum*xsum;
           if (dabs(denom)<MINREAL) {  /* all values at one x point */
              fx[n] = zsum/((double) nsum);
           } else {
              m = (nsum*xzsum - zsum*xsum)/denom;
              b = (zsum - m*xsum)/((double) nsum);

              fx[n] = m*x[n] + b;
              ++nfiltered;
           }
        }
     }
  }

  printf("BestFitLineFilter: %d of %d points are filtered\n",nfiltered,nn);
}

#undef MINREAL

/* ===========================================================================
   BestFitLineDiff (Jun 24/96):

   Routine to take function with finely (irregularly) spaced abscissae, and
   return derivative of function defined with coarse, regularly spaced 
   abscissae.  This function works by finding the best fit line over 
   each interval defined by the coarse grid and returning the slope of
   line.
   ===========================================================================
*/
void BestFitLineDiff(
   double *xraw, double *fxraw, int nraw,
   double lx,double range,
   double *x, double *fx, int nn)
{
  int i,nsum,n;
  double xsum,zsum,xxsum,xzsum,m,hfrng;

  hfrng = 0.5*range;

  for (n=1;n<=nn;n++) {
     nsum = 0;
     zsum = xsum = xzsum = xxsum = 0.0;
     for (i=1;i<=nraw;i++) {
        if (xraw[i]>=x[n]-hfrng && xraw[i]<=x[n]+hfrng) {
           xsum += xraw[i];
           zsum += fxraw[i];
           xzsum += xraw[i]*fxraw[i];
           xxsum += xraw[i]*xraw[i];
           nsum++;
        }
     }

     if (nsum>0) {
        m = (nsum*xzsum - zsum*xsum)/(nsum*xxsum - xsum*xsum);
        fx[n] = m; 
     } else {
        fx[n] = (n==1 ? 0.0 : fx[n-1]);
     }
  }
}
