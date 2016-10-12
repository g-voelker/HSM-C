/* read data from XPlot file */ 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include "myerror.h"
#include "constants.h"
#include "macros.h"

#define MAXCHARLEN 40

/* ======================================================================
   READXPLOT
   ======================================================================
*/
void readXPlot(sfpr, x,fx,ncurves,npts,xzrng, NC,NN)
   int NC,NN;
   int *ncurves,*npts;
   double *x,*fx,*xzrng;
   char *sfpr;
{
   FILE *fpr;
   int i,n,ind,ntotal;
   char filename[MAXFILENAMELEN];
   char dum[MAXCHARLEN];

   printf("Enter XPlot file name: "); fflush(NULL);
   scanf("%s", filename);
   sprintf(sfpr,"%s",filename);
   fpr = fopen(sfpr,"r");


   fscanf(fpr,"%s", dum);
   fscanf(fpr,"%lf",&xzrng[1]);
   fscanf(fpr,"%lf",&xzrng[2]);
   fscanf(fpr,"%lf",&xzrng[3]);
   fscanf(fpr,"%lf",&xzrng[4]);
   fscanf(fpr,"%d", ncurves);

   if (*ncurves>NC) {
      printf( "READXPLOT (ERROR): ncurves=%d too large\n",*ncurves );
      exit(1);
   }
   ntotal=0;
   for (i=1;i<=(*ncurves);i++) {
      fscanf(fpr,"%d", &npts[i] );
      if (ntotal+npts[i]>NN) {
        printf( "READXPLOT (ERROR): npts=%d too large\n",npts[i] );
        exit(1);
      }
      for (n=1;n<=npts[i];n++) {
         ind = ntotal+n;
         fscanf(fpr,"%lf %lf", &x[ind], &fx[ind] );
      }
      ntotal+=npts[i];
   }

   fclose(fpr);
}

/* ======================================================================
   READXPLOT_NOPROMPT: Apr 98

   Read XPlot file without prompting for file name
   ======================================================================
*/
void readXPlot_noprompt(sfpr, x,fx,ncurves,npts,xzrng, NC,NN)
   int NC,NN;
   int *ncurves,*npts;
   double *x,*fx,*xzrng;
   char *sfpr;
{
   FILE *fpr;
   int i,n,ind,ntotal;
   char dum[MAXCHARLEN];

   /* printf( "READXPLOT_NOPROMPT: reading data in file %s\n",sfpr); */

   fpr = fopen(sfpr,"r");

   fscanf(fpr,"%s", dum);
   fscanf(fpr,"%lf",&xzrng[1]);
   fscanf(fpr,"%lf",&xzrng[2]);
   fscanf(fpr,"%lf",&xzrng[3]);
   fscanf(fpr,"%lf",&xzrng[4]);
   fscanf(fpr,"%d", ncurves);

   if (*ncurves>NC) {
      printf( "READXPLOT_NOPROMPT (ERROR): ncurves=%d too large\n",*ncurves );
      exit(1);
   }
   ntotal=0;
   for (i=1;i<=(*ncurves);i++) {
      fscanf(fpr,"%d", &npts[i] );
      for (n=1;n<=npts[i];n++) {
         ind = ntotal+n;
         fscanf(fpr,"%lf %lf", &x[ind], &fx[ind] );
      }
      ntotal+=npts[i];
   }
   if (ntotal>NN) {
      printf( "READXPLOT_NOPROMPT (ERROR): npts=%d too large\n",ntotal);
      exit(1);
   }

   fclose(fpr);
}

/* ===========================================================================
   READXYPLOT: AUG 31/93

   Read data in file save in XYPlot format.
   
   An XYPlot file has the format:

      data.cpt    <filename>
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
void readXYPlot(sfpr, x,y,data,nc,nxypts, NC_XY,NN_XY)
   int NC_XY,NN_XY;
   int *nc;
   int **nxypts;
   double *x,*y,*data;
   char *sfpr;
{
   FILE *fpr;
   int n,i,nx,ny,index,xindex,yindex,ntotal,nptstotal;
   double xmin,xmax,ymin,ymax,dx,dy;
   double clow,chigh,cinc;
   char filename[MAXFILENAMELEN];
   char dum[MAXCHARLEN];

   printf("Enter XYPlot file name: ");
   scanf("%s",filename);
   sprintf(sfpr,"%s",filename); 
   fpr = fopen(sfpr,"r");

   /* ------------------------------------------------------------------
      Read header information.  Ignore values of contour levels.
      ------------------------------------------------------------------
   */
   fscanf(fpr,"%s\n", dum);
   fscanf(fpr,"%lf\n%lf\n%lf\n%lf\n",&xmin,&xmax,&ymin,&ymax);
   fscanf(fpr,"%d", nc);

   if (*nc>NC_XY) {
      printf( "READXYPLOT (ERROR): nc=%d too large\n", *nc );
      exit(1);
   }

   index=xindex=yindex=1;
   for (n=1;n<=(*nc);n++) {
      fscanf(fpr,"%d", &nxypts[n][1] );
      fscanf(fpr,"%d", &nxypts[n][2] );
      fscanf(fpr,"%lf %lf %lf", &clow, &chigh, &cinc);
    
      nx = nxypts[n][1];
      ny = nxypts[n][2];
      ntotal = nx*ny;
      for (i=0;i<ntotal;i++) {
         fscanf(fpr,"%lf", &data[index+i] );
      }

      dx = (xmax-xmin)/((double) (nx-1));
      dy = (ymax-ymin)/((double) (ny-1));
      for (i=0;i<nx;i++) x[xindex+i] = xmin + i*dx;
      for (i=0;i<ny;i++) y[yindex+i] = ymin + i*dy;

      index += ntotal;
      xindex += nx;
      yindex += ny;
   }
   nptstotal=index-1;
   if (nptstotal>NN_XY) {
      printf( "READXYPLOT (ERROR): npts_total=%d too large\n",nptstotal );
      exit(1);
   }

   fclose(fpr);
}

/* ===========================================================================
   READXYPLOT: AUG 31/93

   Read data in file save in XYPlot format.
   
   An XYPlot file has the format:

      data.cpt    <filename>
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
void readXYPlot_noprompt(sfpr, x,y,data,nc,nxypts, NC_XY,NN_XY)
   int NC_XY,NN_XY;
   int *nc;
   int **nxypts;
   double *x,*y,*data;
   char *sfpr;
{
   FILE *fpr;
   int n,i,nx,ny,index,xindex,yindex,ntotal,nptstotal;
   double xmin,xmax,ymin,ymax,dx,dy;
   double clow,chigh,cinc;
   char filename[MAXFILENAMELEN];
   char dum[MAXCHARLEN];

/*
   printf("Enter XYPlot file name: ");
   scanf("%s",filename);
   sprintf(sfpr,"%s",filename); 
*/
   fpr = fopen(sfpr,"r");

   /* ------------------------------------------------------------------
      Read header information.  Ignore values of contour levels.
      ------------------------------------------------------------------
   */
   fscanf(fpr,"%s\n", dum);
   fscanf(fpr,"%lf\n%lf\n%lf\n%lf\n",&xmin,&xmax,&ymin,&ymax);
   fscanf(fpr,"%d", nc);

   if (*nc>NC_XY) {
      printf( "READXYPLOT (ERROR): nc=%d too large\n", *nc );
      exit(1);
   }

   index=xindex=yindex=1;
   for (n=1;n<=(*nc);n++) {
      fscanf(fpr,"%d", &nxypts[n][1] );
      fscanf(fpr,"%d", &nxypts[n][2] );
      fscanf(fpr,"%lf %lf %lf", &clow, &chigh, &cinc);
    
      nx = nxypts[n][1];
      ny = nxypts[n][2];
      ntotal = nx*ny;
      for (i=0;i<ntotal;i++) {
         fscanf(fpr,"%lf", &data[index+i] );
      }

      dx = (xmax-xmin)/((double) (nx-1));
      dy = (ymax-ymin)/((double) (ny-1));
      for (i=0;i<nx;i++) x[xindex+i] = xmin + i*dx;
      for (i=0;i<ny;i++) y[yindex+i] = ymin + i*dy;

      index += ntotal;
      xindex += nx;
      yindex += ny;
   }
   nptstotal=index-1;
   if (nptstotal>NN_XY) {
      printf( "READXYPLOT (ERROR): npts_total=%d too large\n",nptstotal );
      exit(1);
   }

   fclose(fpr);
}

/* ===========================================================================
   READCONPLOT:

   Save data in file for printing by ConPlot routine. 
   This plot format is obsolete but remains for posterity.  See readXYPlot
   for revised version.
   ===========================================================================
*/
void readConPlot(sfpr,fld,nxpts,nzpts, params)
   int *nxpts,*nzpts;
   double **fld,*params;
   char *sfpr;
{
   FILE *fpr;
   int nplt,m,n;
   double alpha,D,ri,re,pr,dnbeta,dmm,dnn,zmin,zmax,xmin,xmax,t;
   double minval,maxval,cinc;
   char dum[MAXFILENAMELEN];

   printf("Enter ConPlot file name:\n");
   scanf("%s",dum);
   sprintf(sfpr,"%s",dum);
   fpr = fopen(sfpr,"r");
   printf("Reading data in file: %s\n", sfpr);

   fscanf(fpr,"%lf\n", &t);
   fscanf(fpr,"%lf\n%lf\n%lf\n", &dmm, &dnn, &dnbeta );
   fscanf(fpr,"%lf\n%lf\n", &D, &alpha);
   fscanf(fpr,"%lf\n%lf\n%lf\n", &ri, &re, &pr );
   fscanf(fpr,"%lf\n%lf\n", &zmin,&zmax );
   fscanf(fpr,"%lf\n%lf\n", &xmin,&xmax );
   fscanf(fpr,"%d\n%d\n", nzpts,nxpts);

   fscanf(fpr,"%d\n", &nplt );
   fscanf(fpr,"%s\n",dum);

   fscanf(fpr,"%lf\n%lf\n%lf\n", &minval,&maxval,&cinc);

   for (m=0;m<= *nzpts-1;m++)
   for (n=0;n<= *nxpts-1;n++) fscanf(fpr,"%lf\n", &fld[n][m] );

   params[2]=alpha;
   params[4]=D;
   params[10]=ri;
   params[11]=re;
   params[12]=pr;
   params[19]= dnbeta;
   params[20]=dmm;
   params[21]=dnn;
   params[22]=zmin;
   params[23]=zmax-zmin;
   params[24]=xmin;
   params[25]=xmax-xmin;
   params[30]=t;
}
