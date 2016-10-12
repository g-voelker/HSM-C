/* ===========================================================================
   Routines to store data in a specified file for plotting by XPlot,XYPlot,
   or ConPlot programs
   ===========================================================================
*/
#define Ncontours 10
#define Minlog -40.0
#define ConMinlog -25.0

#include <stdio.h>
#include <math.h>
#include "complex.h"
#include "alloc_space.h"
#include "constants.h"
#include "macros.h"

/* ===========================================================================
   WRITEXPLOT: AUG 27/93

   XPlot draws a single curve or multiple curves on a graph.  An XPlot file
   has the format:

      data.xpt    <filename>
      xmin
      xmax
      ymin
      ymax
      nc          <number of curves>
      np1         <number of points for curve 1>
      x_1	y_1
      x_2	y_2
      ..    ..
      x_np1	y_np1
      np2         <number of points for curve 2>
      ...etc
     
   ===========================================================================
*/
void writeXPlot(sfpw,x,fnx,ncurves,npts)
   int ncurves;
   int *npts;
   double *x,*fnx;
   char *sfpw;
{
   FILE *fpw;
   int n,i,ind,ntotal;
   double tmp,xmin,xmax,fnmin,fnmax;

   fpw = fopen( sfpw, "w" );

   xmin = xmax = x[1];
   fnmin = fnmax = fnx[1];

   ntotal=0;
   for (i=1;i<=ncurves;i++) {
      for (n=1;n<=npts[i];n++) {
         ind=ntotal+n;
         tmp=x[ind];
         if (xmin>tmp) xmin=tmp;
         if (xmax<tmp) xmax=tmp;
         tmp=fnx[ind];
         if (fnmin>tmp) fnmin=tmp;
         if (fnmax<tmp) fnmax=tmp;
      }
      ntotal+=npts[i];
   }

   fprintf(fpw,"%s\n", sfpw );
   fprintf(fpw,"%15.8e\n%15.8e\n", xmin,xmax );
   fprintf(fpw,"%15.8e\n%15.8e\n", fnmin, fnmax );
   fprintf(fpw,"%d\n", ncurves);

   ntotal=0;
   for (i=1;i<=ncurves;i++) {
      fprintf(fpw,"%d\n", npts[i]);
      for (n=1;n<=npts[i];n++) {
         ind = ntotal+n;
         fprintf(fpw,"%15.8e\t%15.8e\n", x[ind], fnx[ind] );
      }
      ntotal+=npts[i];
   }

   fclose(fpw);
}

/* ===========================================================================
   WRITELOGXPLOT: AUG 27/93

   Write plot in XPlot format but take logarithm of ordinate
   ===========================================================================
*/
void writeLogXPlot(sfpw,x,fnx,ncurves,npts)
   int ncurves;
   int *npts;
   double *x,*fnx;
   char *sfpw;
{
   FILE *fpw;
   int i,n,ntotal,ind;
   double tmp,xmin,xmax,fnmin,fnmax;
   double expMinlog;

   fpw = fopen( sfpw, "w" );

   expMinlog = exp(Minlog);
   if (fnx[1]<=0.0) fnx[1] = expMinlog;
   xmin = xmax = x[1];
   fnmin = fnmax = log(fnx[1]);

   ntotal=0;
   for (i=1;i<=ncurves;i++) {
      for (n=1;n<=npts[i];n++) {
         ind=ntotal+n;
         tmp=x[ind];
         if (xmin>tmp) xmin=tmp;
         if (xmax<tmp) xmax=tmp;
         if (fnx[ind]<=0.0) fnx[ind]=expMinlog;
         tmp=log(fnx[ind]);
         if (fnmin>tmp) fnmin=tmp;
         if (fnmax<tmp) fnmax=tmp;
      }
      ntotal+=npts[i];
   }

   fprintf(fpw,"%s\n", sfpw );
   fprintf(fpw,"%11.4e\n%11.4e\n", xmin,xmax );
   fprintf(fpw,"%11.4e\n%11.4e\n", fnmin, fnmax );
   fprintf(fpw,"%d\n", ncurves);

   ntotal=0;
   for (i=1;i<=ncurves;i++) {
      fprintf(fpw,"%d\n", npts[i]);
      for (n=1;n<=npts[i];n++) { 
         ind=ntotal+n;
         fprintf(fpw,"%11.4e\t%11.4e\n", x[ind], log(fnx[ind]) );
      }
      ntotal+=npts[i];
   }

   fclose(fpw);
}

/* ===========================================================================
   WRITELOGLOGXPLOT: AUG 27/93

   Write plot in XPlot format but take logarithm of abcissae and ordinate
   ===========================================================================
*/
void writeLogLogXPlot(sfpw,x,fnx,ncurves,npts)
   int ncurves;
   int *npts;
   double *x,*fnx;
   char *sfpw;
{
   FILE *fpw;
   int i,n,ntotal,ind;
   double tmp,expMinlog,xmin,xmax,fnmin,fnmax;

   fpw = fopen( sfpw, "w" );

   expMinlog = exp(Minlog);
   if (x[1]<=0.0) x[1] = expMinlog;
   if (fnx[1]<=0.0) fnx[1] = expMinlog;
   xmin = xmax = log(x[1]);
   fnmin = fnmax = log(fnx[1]);

   ntotal=0;
   for (i=1;i<=ncurves;i++) {
      for (n=2;n<=npts[i];n++) {
         ind = ntotal+n;
         if (x[ind]<=0.0) x[ind]=expMinlog;
         tmp=log(x[ind]);
         if (xmin>tmp) xmin=tmp;
         if (xmax<tmp) xmax=tmp;

         if (fnx[ind]<=0.0) fnx[ind]=expMinlog;
         tmp=log(fnx[ind]);
         if (fnmin>tmp) fnmin=tmp;
         if (fnmax<tmp) fnmax=tmp;
      }
      ntotal+=npts[i];
   }

   fprintf(fpw,"%s\n", sfpw );
   fprintf(fpw,"%11.4e\n%11.4e\n", xmin,xmax );
   fprintf(fpw,"%11.4e\n%11.4e\n", fnmin, fnmax );
   fprintf(fpw,"%d\n", ncurves);

   ntotal=0;
   for (i=1;i<=ncurves;i++) {
      fprintf(fpw,"%d\n", npts[i]);
      for (n=1;n<=npts[i];n++) { 
         ind=ntotal+n;
         fprintf(fpw,"%11.4e\t%11.4e\n", log(x[ind]), log(fnx[ind]) );
      }
      ntotal+=npts[i];
   }

   fclose(fpw);
}

/* ===========================================================================
   WRITEXYPLOT: AUG 31/93

   Save data in file for printing by XYPlot routine. 
   
   An XYPlot file has the format:

      data.xyp    <filename>
      xmin        <span and range of plot>
      xmax
      ymin
      ymax
      nc          <number of contour plots>
      nx1         <number of x and y points for 1st contour plot>
      ny1         
      cllow       <contour intervals>
      clhigh
      clint
      f1_x1_y1    <data>
      f1_x2_y1
      ..
      f1_xnx1_y1
      f1_x1_y2
      ..
      nx2         <number of x and y points for 2nd contour plot>
      ny2         
      c2low       <contour intervals>
      c2high
      c2int
      f2_x1_y1
      ...etc

   Note that if cint<0.0 then the following data is taken to represent
   header information for the plot.
   ===========================================================================
*/

void writeXYPlot(sfpw,x,y,data,nc,nxypts,ctype)
   int nc;
   int *ctype;
   int **nxypts;
   double *x,*y,*data;
   char *sfpw;
{
   FILE *fpw;
   int n,i,nx,ny,npts,index,xindex,yindex;
   double xmin,xmax,ymin,ymax;
   double dmin,dmax,d,clow,chigh,cinc;

   /* ------------------------------------------------------------------
      Determine largest range for xmin/xmax and ymin/ymax

      For now, suppose points for _all_ contour plots are evenly spaced 
      between largest min/max
      ------------------------------------------------------------------
   */
   nx = nxypts[1][1];
   ny = nxypts[1][2];
   xmin = x[1];
   xmax = x[nx];
   ymin = y[1];
   ymax = y[ny];

   xindex=1;
   yindex=1;
   for (n=2;n<=nc;n++) {
      nx = nxypts[n][1];
      ny = nxypts[n][2];
      xmin = min(xmin,x[xindex]);
      xmax = max(xmax,x[xindex+nx-1]);
      ymin = min(ymin,y[yindex]);
      ymax = max(ymax,y[yindex+ny-1]);
      xindex += nx;
      yindex += ny;
   }

   fpw = fopen( sfpw, "w" );

   fprintf(fpw,"%s\n",sfpw);

   fprintf(fpw,"%12.5e\n%12.5e\n", xmin,xmax );
   fprintf(fpw,"%12.5e\n%12.5e\n", ymin,ymax );

   fprintf(fpw,"%d\n", nc);

   index=1;
   for (n=1;n<=nc;n++) {

      nx = nxypts[n][1];
      ny = nxypts[n][2];
      fprintf(fpw,"%d\n%d\n", nx,ny);
      
      npts=nx*ny;
      if (ctype[n]==0) {   /* ignore min/max info or use to print plot info */
         clow = chigh = 0.0;
         cinc = -1.0;
      } else if (ctype[n]>0) {
         dmin = dmax = data[index];
         for (i=1;i<npts;i++) { 
            d = data[index+i];
            if (d<dmin) dmin = d;
            if (d>dmax) dmax = d;
         }
         cinc = (dmax-dmin)/((double) ctype[n]);
         clow = dmin-0.5*cinc;
         chigh = dmax+0.5*cinc;
      }

      fprintf(fpw,"%12.5e\n%12.5e\n%12.5e\n", clow,chigh,cinc);

      for (i=0;i<npts;i++) fprintf(fpw,"%13.6e\n", data[index+i]);

      index+=npts;
   }

   fclose(fpw);
}

/* ===========================================================================
   WRITECONPLOT: AUG 27/93

   Save data in file for printing by ConPlot routine.  This 2D plot format 
   is obsolete but is kept for posterity. (See writeXYPlot for recent version).
   
   A ConPlot file has the format:

      t          <spurious information for printing with plot> 
      mm
      nn
      nbeta
      D
      alpha
      ri
      re
      pr
      ymin        <span and range of contour plot>
      ymax
      xmin
      xmax
      nypts       <number of points in span and range>
      nxpts
      nc          <number of curves:  always = 1>
      data.cpt    <filename>
      cntrmin     <contour min/max and interval>
      cntrmax
      cntrintvl
      data_x1_y1
      data_x2_y1
      data_x3_y1
      ..
      data_xnxpts_y1
      data_x1_y2
      .. etc
     
   ===========================================================================
*/
void writeConPlot(sfpw,fld,nxpts,nzpts, params)
   int nxpts,nzpts;
   double *params;
   double **fld;
   char *sfpw;
{
   FILE *fpw;
   int mm,nn,m,n,nbeta;
   double t,D,alpha,ri,re,pr,H,L,zmin,zmax,xmin,xmax;
   double tmp,minval,maxval,cinc;

   alpha = params[2];
   D = params[4];
   ri = params[10];
   re = params[11];
   pr = params[12];
   nbeta = (int) params[19];   
   mm = (int) params[20];
   nn = (int) params[21];
   zmin = params[22];
   H = params[23];
   xmin = params[24]; 
   L = params[25];
   t = params[30];
   zmax = zmin+H;
   xmax = xmin+L;

   minval = maxval = fld[0][0];
   for (m=0;m<nzpts;m++) for (n=0;n<nxpts;n++) {
      if ((tmp=fld[n][m])<minval ) minval = tmp;
      else if (tmp>maxval) maxval = tmp;
   }
   cinc = (maxval-minval)/((double) Ncontours);

   fpw = fopen( sfpw, "w" );
   fprintf(fpw,"%9.6f\n", t);
   fprintf(fpw,"%5.1f\n%5.1f\n%5.1f\n", 
           ((double) mm), ((double) nn), ((double) nbeta) );
   fprintf(fpw,"%5.2f\n%5.2f\n", D, alpha); 
   fprintf(fpw,"%5.3f\n%5.1f\n%5.1f\n", ri, re, pr ); 
   fprintf(fpw,"%9.6f\n%9.6f\n", zmin,zmax );
   fprintf(fpw,"%9.6f\n%9.6f\n", xmin,xmax );
   fprintf(fpw,"%d\n%d\n", nzpts,nxpts);

   fprintf(fpw,"1\n" );
   fprintf(fpw,"%s\n",sfpw);

   fprintf(fpw,"%9.2e\n%9.2e\n%9.2e\n", minval,maxval,cinc);

   for (m=0;m<nzpts;m++) 
   for (n=0;n<nxpts;n++) fprintf(fpw,"%11.4e\n", fld[n][m] );

   fclose(fpw);
}

/* ===========================================================================
   WRITECONPLOT: AUG 27/93

   Save logarithm of data in file for printing by ConPlot routine.
   ===========================================================================
*/
void writeLogConPlot(sfpw,fld,nxpts,nzpts, params)
   int nxpts,nzpts;
   double *params;
   double **fld;
   char *sfpw;
{
   FILE *fpw;
   int mm,nn,m,n,nbeta;
   double t,D,alpha,ri,re,pr,H,L,zmin,zmax,xmin,xmax;
   double tmp,minval,maxval,cinc;
   double expMinlog;

   alpha = params[2];
   D = params[4];
   ri = params[10];
   re = params[11];
   pr = params[12];
   nbeta = (int) params[19];   
   mm = (int) params[20];
   nn = (int) params[21];
   zmin = params[22];
   H = params[23];
   xmin = params[24]; 
   L = params[25];
   t = params[30];
   zmax = zmin+H;
   xmax = xmin+L;

   expMinlog = exp(ConMinlog);

   for (n=0;n<nxpts;n++) 
   for (m=0;m<nzpts;m++) 
      if (dabs(fld[n][m])<expMinlog) fld[n][m] = expMinlog;

   minval = maxval = log(dabs(fld[0][0]));
   for (m=0;m<=nzpts;m++) for (n=0;n<=nxpts;n++) {
      if ((tmp=log(dabs(fld[n][m])))<minval ) minval = tmp;
      else if (tmp>maxval) maxval = tmp;
   }
   cinc = (maxval-minval)/((double) Ncontours);

   fpw = fopen( sfpw, "w" );
   fprintf(fpw,"%9.6f\n", t);
   fprintf(fpw,"%5.1f\n%5.1f\n%5.1f\n", 
           ((double) mm), ((double) nn), ((double) nbeta) );
   fprintf(fpw,"%5.2f\n%5.2f\n", D, alpha); 
   fprintf(fpw,"%5.3f\n%5.1f\n%5.1f\n", ri, re, pr ); 
   fprintf(fpw,"%9.6f\n%9.6f\n", zmin,zmax );
   fprintf(fpw,"%9.6f\n%9.6f\n", xmin,xmax );
   fprintf(fpw,"%d\n%d\n", nzpts,nxpts);

   fprintf(fpw,"1\n" );
   fprintf(fpw,"%s\n",sfpw);

   fprintf(fpw,"%9.2e\n%9.2e\n%9.2e\n", minval,maxval,cinc);

   for (m=0;m<nzpts;m++) 
   for (n=0;n<nxpts;n++) fprintf(fpw,"%9.2e\n", log(dabs(fld[n][m])) );

   fclose(fpw);
}
