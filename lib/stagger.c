/* ============================================================================
   Files to take a field and define intermediate points by linear interpolation
   ============================================================================
*/
#include "complex.h"
#include "alloc_space.h"
#include "constants.h"
#include "macros.h"

/* ===========================================================================
   RSTAGGER: May 3/94

   Take a field defined at evenly spaced points i=0..n
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required
   ===========================================================================
*/
void Rstagger(f,ox,nx, staggerflg, sf)
   int ox,nx,staggerflg;
   double *f,*sf;
{
   int istflg,i;
   int ip1,im1;

   istflg = iabs(staggerflg);

   if (istflg!=1) myerror("RSTAGGER: staggerflg out of bounds");
   if (staggerflg==1 && ox!=1) myerror("RSTAGGER: ox should equal one");
   if (staggerflg==-1 && ox!=0) myerror("RSTAGGER: ox should equal zero");

   if (staggerflg==1) {
      for (i=1;i<nx;i++) {
         ip1 = i+1;
         sf[i] = 0.5*(f[ip1]+f[i]);
      }
      sf[0] = 1.5*f[1] - 0.5*f[2];
      sf[nx] = 1.5*f[nx] - 0.5*f[nx-1];
   } else if (staggerflg==-1) {
      for (i=nx;i>=1;i--) {
         im1 = i-1;
         sf[i] = 0.5*(f[im1]+f[i]);
      }
   }
}

/* ===========================================================================
   _RSTAGGER: May 3/94

   Take a field defined at evenly spaced points i=0..n
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required
   ===========================================================================
*/
void _Rstagger(f,ox,nx, staggerflg)
   int ox,nx,staggerflg;
   double *f;
{
   int istflg,i;
   int ip1,im1;
   double tmp;

   istflg = iabs(staggerflg);

   if (istflg!=1) myerror("RSTAGGER: staggerflg out of bounds");
   if (staggerflg==1 && ox!=1) myerror("RSTAGGER: ox should equal one");
   if (staggerflg==-1 && ox!=0) myerror("RSTAGGER: ox should equal zero");

   if (staggerflg==1) {
      tmp = 1.5*f[nx] - 0.5*f[nx-1];
      for (i=1;i<nx;i++) {
         ip1 = i+1;
         f[i] = 0.5*(f[ip1]+f[i]);
      }
      f[0] = 1.5*f[1] - 0.5*f[2];
      f[nx] = tmp;
   } else if (staggerflg==-1) {
      for (i=nx;i>=1;i--) {
         im1 = i-1;
         f[i] = 0.5*(f[im1]+f[i]);
      }
   }
}

/* ===========================================================================
   RSTAGGER2D: Sept 30/93

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  
   If staggerflg == +/- 1  average in y direction alone.
   If staggerflg == +/- 2  average in z direction alone.

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required

   ===========================================================================
*/
void Rstagger2D(f,oy,ny,oz,nz, staggerflg, sf)
   int oy,ny,oz,nz,staggerflg;
   double **f,**sf;
{
   int istflg,j,k;
   int jp1,jm1,kp1,km1;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>2) 
      myerror("RSTAGGER2D: staggerflg out of bounds");

   if (istflg==1 && ny==0) myerror("RSTAGGER2D: ny equals zero");
   if (istflg==2 && nz==0) myerror("RSTAGGER2D: nz equals zero");
   if (staggerflg==1 && oy!=1) myerror("RSTAGGER2D: oy should equal one");
   if (staggerflg==2 && oz!=1) myerror("RSTAGGER2D: oz should equal one");
   if (staggerflg==-1 && oy!=0) myerror("RSTAGGER2D: oy should equal zero");
   if (staggerflg==-2 && oz!=0) myerror("RSTAGGER2D: oz should equal zero");

   if (staggerflg==1) {
      for (j=1;j<ny;j++) {
         jp1 = j+1;
         for (k=oz;k<=nz;k++) 
            sf[j][k] = 0.5*(f[jp1][k]+f[j][k]);
      }
      for (k=oz;k<=nz;k++) {
         sf[ny][k] = 1.5*f[ny][k] - 0.5*f[ny-1][k];
         sf[0][k] = 1.5*f[1][k] - 0.5*f[2][k];
      }
   } else if (staggerflg==2) {
      for (k=1;k<nz;k++) {
         kp1=k+1;
         for (j=oy;j<=ny;j++) 
            sf[j][k] = 0.5*(f[j][kp1]+f[j][k]);
      }
      for (j=oy;j<=ny;j++) { 
         sf[j][nz] = 1.5*f[j][nz] - 0.5*f[j][nz-1];
         sf[j][0] = 1.5*f[j][1] - 0.5*f[j][2];
      }
      
   } else if (staggerflg==-1) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (k=oz;k<=nz;k++) 
            sf[j][k] = 0.5*(f[jm1][k]+f[j][k]);
      }
      
   } else if (staggerflg==-2) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (j=oy;j<=ny;j++)
            sf[j][k] = 0.5*(f[j][km1]+f[j][k]);
      }
   }
}

/* ===========================================================================
   _RSTAGGER2D: Sept 30/93

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  

   If staggerflg==1  average in y direction alone.
   If staggerflg==2  average in z direction alone.

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required
   ===========================================================================
*/
void _Rstagger2D(f,oy,ny,oz,nz, staggerflg)
   int oy,ny,oz,nz,staggerflg;
   double **f;
{
   int istflg,j,k;
   int jp1,jm1,kp1,km1;
   double tmp;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>2) 
      myerror("_RSTAGGER2D: staggerflg out of bounds");

   if (istflg==1 && ny==0) myerror("_RSTAGGER2D: ny equals zero");
   if (istflg==2 && nz==0) myerror("_RSTAGGER2D: nz equals zero");
   if (staggerflg==1 && oy!=1) myerror("_RSTAGGER2D: oy should equal one");
   if (staggerflg==2 && oz!=1) myerror("_RSTAGGER2D: oz should equal one");
   if (staggerflg==-1 && oy!=0) myerror("_RSTAGGER2D: oy should equal zero");
   if (staggerflg==-2 && oz!=0) myerror("_RSTAGGER2D: oz should equal zero");

   if (staggerflg==1) {
      for (k=oz;k<=nz;k++) { 
         tmp = 1.5*f[ny][k] - 0.5*f[ny-1][k];
         f[0][k] = 1.5*f[1][k] - 0.5*f[2][k];
         for (j=1;j<ny;j++) f[j][k] = 0.5*(f[j+1][k]+f[j][k]);
         f[ny][k] = tmp;
      }
   } else if (staggerflg==2) {
      for (j=oy;j<=ny;j++) {
         tmp = 1.5*f[j][nz] - 0.5*f[j][nz-1];
         f[j][0] = 1.5*f[j][1] - 0.5*f[j][2];
         for (k=1;k<nz;k++) f[j][k] = 0.5*(f[j][k+1]+f[j][k]);
         f[j][nz] = tmp;
      }
   } else if (staggerflg==-1) {
      for (j=ny;j>=1;j--) {
         jm1=j-1; 
         for (k=oz;k<=nz;k++)  f[j][k] = 0.5*(f[jm1][k]+f[j][k]);
      }
   } else if (staggerflg==-2) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (j=oy;j<=ny;j++) f[j][k] = 0.5*(f[j][km1]+f[j][k]);
      }
   }
}

/* ===========================================================================
   RPSTAGGER2D: Jan 14/94

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  
   If staggerflg == +/- 1  average in y direction alone.
   If staggerflg == +/- 2  average in z direction alone.

   Staggering assumes grid range is (0,1,2,..,n) and that f[0],f[n] are
   identical (i.e. is periodic in that direction).

   if staggerflg is positive then x_(i), x_(i+1)  determine x'_(i)
   if staggerflg is negative then x_(i-1), x_(i)  determine x'_(i)

   note  x_0 = x_nx  =>  x_(-1)==x_(nx-1), x_(nx+1)==x_(1)
   ===========================================================================
*/
void RPstagger2D(f,oy,ny,oz,nz, staggerflg,sf)
   int oy,ny,oz,nz,staggerflg;
   double **f,**sf;
{
   int istflg,j,k;
   int jp1,jm1,kp1,km1;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>2) 
      myerror("RPSTAGGER2D: staggerflg out of bounds");

   if (istflg==1 && ny==0) myerror("RPSTAGGER2D: ny equals zero");
   if (istflg==2 && nz==0) myerror("RPSTAGGER2D: nz equals zero");
   if (istflg==1 && oy!=0) myerror("RPSTAGGER2D: oy should equal zero");
   if (istflg==2 && oz!=0) myerror("RPSTAGGER2D: oz should equal zero");

   if (staggerflg==1) {
      for (j=0;j<ny;j++) {
         jp1 = j+1;
         for (k=oz;k<=nz;k++) 
            sf[j][k] = 0.5*(f[jp1][k]+f[j][k]);
      }
      for (k=oz;k<=nz;k++) sf[ny][k] = sf[0][k];
   } else if (staggerflg==2) {
      for (k=0;k<nz;k++) {
         kp1=k+1;
         for (j=oy;j<=ny;j++) 
            sf[j][k] = 0.5*(f[j][kp1]+f[j][k]);
      }
      for (j=oy;j<=ny;j++) sf[j][nz] = sf[j][0];
   } else if (staggerflg==-1) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (k=oz;k<=nz;k++) 
            sf[j][k] = 0.5*(f[jm1][k]+f[j][k]);
      }
      for (k=oz;k<=nz;k++) sf[0][k] = sf[ny][k];
   } else if (staggerflg==-2) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (j=oy;j<=ny;j++)
            sf[j][k] = 0.5*(f[j][km1]+f[j][k]);
      }
      for (j=oy;j<=ny;j++) sf[j][0] = sf[j][nz];
   }
}

/* ===========================================================================
   _RPSTAGGER2D: Jan 14/94

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  

   If staggerflg==1  average in y direction alone.
   If staggerflg==2  average in z direction alone.

   Assume periodic boundary conditions in direction of staggering.
   The stagger index should range from (0,1,2,...n) where 0 and n are
   identified (See RPSTAGGER2D).
   ===========================================================================
*/
void _RPstagger2D(f,oy,ny,oz,nz, staggerflg)
   int oy,ny,oz,nz,staggerflg;
   double **f;
{
   int istflg,j,k;
   int jp1,jm1,kp1,km1;
   double tmp;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>2) 
      myerror("_RPSTAGGER2D: staggerflg out of bounds");

   if (istflg==1 && ny==0) myerror("_RPSTAGGER2D: ny equals zero");
   if (istflg==2 && nz==0) myerror("_RPSTAGGER2D: nz equals zero");
   if (istflg==1 && oy!=0) myerror("_RPSTAGGER2D: oy should equal zero");
   if (istflg==2 && oz!=0) myerror("_RPSTAGGER2D: oz should equal zero");

   if (staggerflg==1) {
      for (j=0;j<ny;j++) {
         jp1 = j+1;
         for (k=oz;k<=nz;k++) f[j][k] = 0.5*(f[jp1][k]+f[j][k]);
      }
      for (k=oz;k<=nz;k++) f[ny][k] = f[0][k];
   } else if (staggerflg==2) {
      for (k=0;k<nz;k++) {
         kp1=k+1;
         for (j=oy;j<=ny;j++) f[j][k] = 0.5*(f[j][kp1]+f[j][k]);
      }
      for (j=oy;j<=ny;j++) f[j][nz] = f[j][0];
   } else if (staggerflg==-1) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (k=oz;k<=nz;k++) f[j][k] = 0.5*(f[jm1][k]+f[j][k]);
      }
      for (k=oz;k<=nz;k++) f[0][k] = f[ny][k];
   } else if (staggerflg==-2) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (j=oy;j<=ny;j++) f[j][k] = 0.5*(f[j][km1]+f[j][k]);
      }
      for (j=oy;j<=ny;j++) f[j][0] = f[j][nz];
   }
}

/* ===========================================================================
   RZSTAGGER2D: Feb 21/94

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  
   If staggerflg == +/- 1  average in y direction alone.
   If staggerflg == +/- 2  average in z direction alone.

   Staggering assumes grid range is (0,1,2,..,n) and that f[0],f[n] are
   identical (i.e. is periodic in that direction).

   if staggerflg is positive then x_(i), x_(i+1)  determine x'_(i)
   if staggerflg is negative then x_(i-1), x_(i)  determine x'_(i)

   note  y_0 = y_ny = 0
   ===========================================================================
*/
void RZstagger2D(f,oy,ny,oz,nz, staggerflg,sf)
   int oy,ny,oz,nz,staggerflg;
   double **f,**sf;
{
   int istflg,j,k;
   int jp1,jm1,kp1,km1;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>2) 
      myerror("RZSTAGGER2D: staggerflg out of bounds");

   if (istflg==1 && ny==0) myerror("RZSTAGGER2D: ny equals zero");
   if (istflg==2 && nz==0) myerror("RZSTAGGER2D: nz equals zero");
   if (staggerflg==1 && oy!=1) myerror("RZSTAGGER2D: oy should equal one");
   if (staggerflg==2 && oz!=1) myerror("RZSTAGGER2D: oz should equal one");
   if (staggerflg==-1 && oy!=0) myerror("RZSTAGGER2D: oy should equal zero");
   if (staggerflg==-2 && oz!=0) myerror("RZSTAGGER2D: oz should equal zero");

   if (staggerflg==1) {
      for (j=1;j<ny;j++) {
         jp1 = j+1;
         for (k=oz;k<=nz;k++) 
            sf[j][k] = 0.5*(f[jp1][k]+f[j][k]);
      }
      for (k=oz;k<=nz;k++) sf[0][k] = sf[ny][k] = 0.0;
   } else if (staggerflg==2) {
      for (k=1;k<nz;k++) {
         kp1=k+1;
         for (j=oy;j<=ny;j++) 
            sf[j][k] = 0.5*(f[j][kp1]+f[j][k]);
      }
      for (j=oy;j<=ny;j++) sf[j][0]=sf[j][nz]=0.0;
   } else if (staggerflg==-1) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (k=oz;k<=nz;k++) 
            sf[j][k] = 0.5*(f[jm1][k]+f[j][k]);
      }
   } else if (staggerflg==-2) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (j=oy;j<=ny;j++)
            sf[j][k] = 0.5*(f[j][km1]+f[j][k]);
      }
   }
}

/* ===========================================================================
   _RZSTAGGER2D: Feb 21/94

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  

   If staggerflg==1  average in y direction alone.
   If staggerflg==2  average in z direction alone.

   Assume zero boundary conditions in direction of staggering.
   The stagger index should range from (0,1,2,...n) where 0 and n are
   identified (See RZSTAGGER2D).
   ===========================================================================
*/
void _RZstagger2D(f,oy,ny,oz,nz, staggerflg)
   int oy,ny,oz,nz,staggerflg;
   double **f;
{
   int istflg,j,k;
   int jp1,jm1,kp1,km1;
   double tmp;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>2) 
      myerror("_RZSTAGGER2D: staggerflg out of bounds");

   if (istflg==1 && ny==0) myerror("_RZSTAGGER2D: ny equals zero");
   if (istflg==2 && nz==0) myerror("_RZSTAGGER2D: nz equals zero");
   if (staggerflg==1 && oy!=1) myerror("_RZSTAGGER2D: oy should equal one");
   if (staggerflg==2 && oz!=1) myerror("_RZSTAGGER2D: oz should equal one");
   if (staggerflg==-1 && oy!=0) myerror("_RZSTAGGER2D: oy should equal zero");
   if (staggerflg==-2 && oz!=0) myerror("_RZSTAGGER2D: oz should equal zero");

   if (staggerflg==1) {
      for (j=1;j<ny;j++) {
         jp1 = j+1;
         for (k=oz;k<=nz;k++) f[j][k] = 0.5*(f[jp1][k]+f[j][k]);
      }
      for (k=oz;k<=nz;k++) f[0][k]=f[ny][k]=0.0;
   } else if (staggerflg==2) {
      for (k=1;k<nz;k++) {
         kp1=k+1;
         for (j=oy;j<=ny;j++) f[j][k] = 0.5*(f[j][kp1]+f[j][k]);
      }
      for (j=oy;j<=ny;j++) f[j][nz]=f[j][0]=0.0;
   } else if (staggerflg==-1) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (k=oz;k<=nz;k++) f[j][k] = 0.5*(f[jm1][k]+f[j][k]);
      }
   } else if (staggerflg==-2) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (j=oy;j<=ny;j++) f[j][k] = 0.5*(f[j][km1]+f[j][k]);
      }
   }
}

/* ===========================================================================
   CSTAGGER: Oct 28/93

   Take a field defined at evenly spaced points i=0..n
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required
   ===========================================================================
*/
void Cstagger(f,ox,nx, staggerflg,sf)
   int ox,nx,staggerflg;
   dcomplex *f,*sf;
{
   int istflg,i;
   int ip1,im1;

   istflg = iabs(staggerflg);

   if (istflg!=1) myerror("CSTAGGER: staggerflg out of bounds");
   if (staggerflg==1 && ox!=1) myerror("CSTAGGER: ox should equal one");
   if (staggerflg==-1 && ox!=0) myerror("CSTAGGER: ox should equal zero");

   if (staggerflg==1) {
      for (i=1;i<nx;i++) {
         ip1 = i+1;
         sf[i].r = 0.5*(f[ip1].r+f[i].r);
         sf[i].i = 0.5*(f[ip1].i+f[i].i);
      }
      sf[0].r = 1.5*f[1].r - 0.5*f[2].r;
      sf[0].i = 1.5*f[1].i - 0.5*f[2].i;
      sf[nx].r = 1.5*f[nx].r - 0.5*f[nx-1].r;
      sf[nx].i = 1.5*f[nx].i - 0.5*f[nx-1].i;
   } else if (staggerflg==-1) {
      for (i=nx;i>=1;i--) {
         im1 = i-1;
         sf[i].r = 0.5*(f[im1].r+f[i].r);
         sf[i].i = 0.5*(f[im1].i+f[i].i);
      }
   }
}

/* ===========================================================================
   CSTAGGER2D: Sept 30/93

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  
   If staggerflg==1  average in y direction alone.
   If staggerflg==2  average in z direction alone.

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required
   ===========================================================================
*/
void Cstagger2D(f,oy,ny,oz,nz, staggerflg,sf)
   int oy,ny,oz,nz,staggerflg;
   dcomplex **f,**sf;
{
   int istflg,j,k;
   int jp1,jm1,kp1,km1;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>2) 
      myerror("CSTAGGER2D: staggerflg out of bounds");

   if (istflg==1 && ny==0) myerror("CSTAGGER2D: ny equals zero");
   if (istflg==2 && nz==0) myerror("CSTAGGER2D: nz equals zero");
   if (staggerflg==1 && oy!=1) myerror("CSTAGGER2D: oy should equal one");
   if (staggerflg==2 && oz!=1) myerror("CSTAGGER2D: oz should equal one");
   if (staggerflg==-1 && oy!=0) myerror("CSTAGGER2D: oy should equal zero");
   if (staggerflg==-2 && oz!=0) myerror("CSTAGGER2D: oz should equal zero");

   if (staggerflg==1) {
      for (j=1;j<ny;j++) {
         jp1 = j+1;
         for (k=oz;k<=nz;k++) {
            sf[j][k].r = 0.5*(f[jp1][k].r+f[j][k].r);
            sf[j][k].i = 0.5*(f[jp1][k].i+f[j][k].i);
         }
      }
      for (k=oz;k<=nz;k++) {
         sf[0][k].r = 1.5*f[1][k].r - 0.5*f[2][k].r;
         sf[0][k].i = 1.5*f[1][k].i - 0.5*f[2][k].i;
         sf[ny][k].r = 1.5*f[ny][k].r - 0.5*f[ny-1][k].r;
         sf[ny][k].i = 1.5*f[ny][k].i - 0.5*f[ny-1][k].i;
      }
   } else if (staggerflg==2) {
      for (k=1;k<nz;k++) {
         kp1=k+1;
         for (j=oy;j<=ny;j++) {
            sf[j][k].r = 0.5*(f[j][kp1].r+f[j][k].r);
            sf[j][k].i = 0.5*(f[j][kp1].i+f[j][k].i);
         }
      }
      for (j=oy;j<=ny;j++) {
         sf[j][0].r = 1.5*f[j][1].r - 0.5*f[j][2].r;
         sf[j][0].i = 1.5*f[j][1].i - 0.5*f[j][2].i;
         sf[j][nz].r = 1.5*f[j][nz].r - 0.5*f[j][nz-1].r;
         sf[j][nz].i = 1.5*f[j][nz].i - 0.5*f[j][nz-1].i;
      }
   } else if (staggerflg==-1) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (k=oz;k<=nz;k++) {
            sf[j][k].r = 0.5*(f[jm1][k].r+f[j][k].r);
            sf[j][k].i = 0.5*(f[jm1][k].i+f[j][k].i);
         }
      }
   } else if (staggerflg==-2) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (j=oy;j<=ny;j++) {
            sf[j][k].r = 0.5*(f[j][km1].r+f[j][k].r);
            sf[j][k].i = 0.5*(f[j][km1].i+f[j][k].i);
         }
      }
   }
}

/* ===========================================================================
   _CSTAGGER2D: Sept 30/93

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  
   If staggerflg==1  average in y direction alone.
   If staggerflg==2  average in z direction alone.

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required
   ===========================================================================
*/
void _Cstagger2D(f,oy,ny,oz,nz, staggerflg)
   int oy,ny,oz,nz,staggerflg;
   dcomplex **f;
{
   int istflg,j,k;
   int jp1,jm1,kp1,km1;
   dcomplex ctmp;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>2) 
      myerror("_CSTAGGER: staggerflg out of bounds");

   if (istflg==1 && ny==0) myerror("_CSTAGGER2D: ny equals zero");
   if (istflg==2 && nz==0) myerror("_CSTAGGER2D: nz equals zero");
   if (staggerflg==1 && oy!=1) myerror("CSTAGGER2D: oy should equal one");
   if (staggerflg==2 && oz!=1) myerror("CSTAGGER2D: oz should equal one");
   if (staggerflg==-1 && oy!=0) myerror("CSTAGGER2D: oy should equal zero");
   if (staggerflg==-2 && oz!=0) myerror("CSTAGGER2D: oz should equal zero");

   if (staggerflg==1) {
      for (k=oz;k<=nz;k++) {
         ctmp.r = 1.5*f[ny][k].r - 0.5*f[ny-1][k].r;
         ctmp.i = 1.5*f[ny][k].i - 0.5*f[ny-1][k].i;
         f[0][k].r = 1.5*f[1][k].r - 0.5*f[2][k].r;
         f[0][k].i = 1.5*f[1][k].i - 0.5*f[2][k].i;
         for (j=1;j<ny;j++) {
            jp1 = j+1;
            f[j][k].r = 0.5*(f[jp1][k].r+f[j][k].r);
            f[j][k].i = 0.5*(f[jp1][k].i+f[j][k].i);
         }
         f[ny][k] = ctmp;
      }
   } else if (staggerflg==2) {
      for (j=oy;j<=ny;j++) {
         ctmp.r = 1.5*f[j][nz].r - 0.5*f[j][nz-1].r;
         ctmp.i = 1.5*f[j][nz].i - 0.5*f[j][nz-1].i;
         f[j][0].r = 1.5*f[j][1].r - 0.5*f[j][2].r;
         f[j][0].i = 1.5*f[j][1].i - 0.5*f[j][2].i;
         for (k=1;k<nz;k++) {
            kp1 = k+1;
            f[j][k].r = 0.5*(f[j][kp1].r+f[j][k].r);
            f[j][k].i = 0.5*(f[j][kp1].i+f[j][k].i);
         }
         f[j][nz] = ctmp;
      }
   } else if (staggerflg==-1) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (k=oz;k<=nz;k++) {
            f[j][k].r = 0.5*(f[jm1][k].r+f[j][k].r);
            f[j][k].i = 0.5*(f[jm1][k].i+f[j][k].i);
         }
      }
   } else if (staggerflg==-2) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (j=oy;j<=ny;j++) {
            f[j][k].r = 0.5*(f[j][km1].r+f[j][k].r);
            f[j][k].i = 0.5*(f[j][km1].i+f[j][k].i);
         }
      }
   }
}

/* ===========================================================================
   RSTAGGER3D: Nov 25/93

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  
   If staggerflg==1  average in x and z directions.
   If staggerflg==2  average in y direction alone.
   If staggerflg==3  average in z direction alone.

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required
   ===========================================================================
*/
void Rstagger3D(f,ox,nx,oy,ny,oz,nz, staggerflg,sf)
   int ox,nx,oy,ny,oz,nz,staggerflg;
   double ***f,***sf;
{
   int istflg,n,j,k;
   int np1,nm1,jp1,jm1,kp1,km1;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>3) 
      myerror("RSTAGGER3D: staggerflg out of bounds");

   if (staggerflg==1 && ox!=1) myerror("RSTAGGER3D: ox should equal one");
   if (staggerflg==2 && oy!=1) myerror("RSTAGGER3D: oy should equal one");
   if (staggerflg==3 && oz!=1) myerror("RSTAGGER3D: oz should equal one");
   if (staggerflg==-1 && ox!=0) myerror("RSTAGGER3D: ox should equal zero");
   if (staggerflg==-2 && oy!=0) myerror("RSTAGGER3D: oy should equal zero");
   if (staggerflg==-3 && oz!=0) myerror("RSTAGGER3D: oz should equal zero");

   if (staggerflg==1) {
      for (n=1;n<nx;n++) {
         np1=n+1;
         for (j=oy;j<=ny;j++) 
         for (k=oz;k<=nz;k++) 
            sf[n][j][k] = 0.5*(f[np1][j][k]+f[n][j][k]);
      }
      for (j=oy;j<=ny;j++) 
      for (k=oz;k<=nz;k++) {
         sf[0][j][k] = 1.5*f[1][j][k] - 0.5*f[2][j][k];
         sf[nx][j][k] = 1.5*f[nx][j][k] - 0.5*f[nx-1][j][k];
      }

   } else if (staggerflg==2) {
      for (j=1;j<ny;j++) {
         jp1 = j+1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) 
            sf[n][j][k] = 0.5*(f[n][jp1][k]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (k=oz;k<=nz;k++) { 
         sf[n][0][k] = 1.5*f[n][1][k] - 0.5*f[n][2][k];
         sf[n][ny][k] = 1.5*f[n][ny][k] - 0.5*f[n][ny-1][k];
      }
   } else if (staggerflg==3) {
      for (k=1;k<nz;k++) {
         kp1=k+1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) 
            sf[n][j][k] = 0.5*(f[n][j][kp1]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (j=oy;j<=ny;j++) { 
         sf[n][j][0] = 1.5*f[n][j][1] - 0.5*f[n][j][2];
         sf[n][j][nz] = 1.5*f[n][j][nz] - 0.5*f[n][j][nz-1];
      }
   } else if (staggerflg==-1) {
      for (n=nx;n>=1;n--) {
         nm1=n-1;
         for (j=oy;j<=ny;j++) 
         for (k=oz;k<=nz;k++) 
            sf[n][j][k] = 0.5*(f[nm1][j][k]+f[n][j][k]);
      }
   } else if (staggerflg==-2) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) 
            sf[n][j][k] = 0.5*(f[n][jm1][k]+f[n][j][k]);
      }
   } else if (staggerflg==-3) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) 
            sf[n][j][k] = 0.5*(f[n][j][km1]+f[n][j][k]);
      }
   }
}

/* ===========================================================================
   _RSTAGGER3D: Nov 25/93

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  
   If staggerflg==1  average in x direction alone.
   If staggerflg==2  average in y direction alone.
   If staggerflg==3  average in z direction alone.

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required
   ===========================================================================
*/
void _Rstagger3D(f,ox,nx,oy,ny,oz,nz, staggerflg)
   int ox,nx,oy,ny,oz,nz,staggerflg;
   double ***f;
{
   int istflg,n,j,k;
   int np1,nm1,jp1,jm1,kp1,km1;
   double tmp;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>3) 
      myerror("_RSTAGGER3D: staggerflg out of bounds");

   if (staggerflg==1 && ox!=1) myerror("_RSTAGGER3D: ox should equal one");
   if (staggerflg==2 && oy!=1) myerror("_RSTAGGER3D: oy should equal one");
   if (staggerflg==3 && oz!=1) myerror("_RSTAGGER3D: oz should equal one");
   if (staggerflg==-1 && ox!=0) myerror("_RSTAGGER3D: ox should equal zero");
   if (staggerflg==-2 && oy!=0) myerror("_RSTAGGER3D: oy should equal zero");
   if (staggerflg==-3 && oz!=0) myerror("_RSTAGGER3D: oz should equal zero");

   if (staggerflg==1) {
      for (j=oy;j<=ny;j++) 
      for (k=oz;k<=nz;k++) { 
         tmp = 1.5*f[nx][j][k] - 0.5*f[nx-1][j][k];
         f[0][j][k] = 1.5*f[1][j][k] - 0.5*f[2][j][k];
         for (n=1;n<nx;n++) f[n][j][k] = 0.5*(f[n+1][j][k]+f[n][j][k]);
         f[nx][j][k] = tmp;
      }
   } else if (staggerflg==2) {
      for (n=ox;n<=nx;n++)
      for (k=oz;k<=nz;k++) { 
         tmp = 1.5*f[n][ny][k] - 0.5*f[n][ny-1][k];
         f[n][0][k] = 1.5*f[n][1][k] - 0.5*f[n][2][k];
         for (j=1;j<ny;j++) f[n][j][k] = 0.5*(f[n][j+1][k]+f[n][j][k]);
         f[n][ny][k] = tmp;
      }
   } else if (staggerflg==3) {
      for (n=ox;n<=nx;n++)
      for (j=oy;j<=ny;j++) { 
         tmp = 1.5*f[n][j][nz] - 0.5*f[n][j][nz-1];
         f[n][j][0] = 1.5*f[n][j][1] - 0.5*f[n][j][2];
         for (k=1;k<nz;k++) f[n][j][k] = 0.5*(f[n][j][k+1]+f[n][j][k]);
         f[n][j][nz] = tmp;
      }
   } else if (staggerflg==-1) {
      for (n=nx;n>=1;n--) {
         nm1=n-1;
         for (j=oy;j<=ny;j++) 
         for (k=oz;k<=nz;k++) 
            f[n][j][k] = 0.5*(f[nm1][j][k]+f[n][j][k]);
      }
   } else if (staggerflg==-2) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) 
            f[n][j][k] = 0.5*(f[n][jm1][k]+f[n][j][k]);
      }
   } else if (staggerflg==-3) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) 
            f[n][j][k] = 0.5*(f[n][j][km1]+f[n][j][k]);
      }
   }
}

/* ===========================================================================
   RPSTAGGER3D: Jan 14/94

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  
   If staggerflg==1  average in x and z directions.
   If staggerflg==2  average in y direction alone.
   If staggerflg==3  average in z direction alone.

   This routine assumes periodicity in the staggerflg direction with
   index ranging from (0,1,2,..,n).  (See RPSTAGGER2D)
   ===========================================================================
*/
void RPstagger3D(f,ox,nx,oy,ny,oz,nz, staggerflg,sf)
   int ox,nx,oy,ny,oz,nz,staggerflg;
   double ***f,***sf;
{
   int istflg,n,j,k;
   int np1,nm1,jp1,jm1,kp1,km1;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>3) myerror("RPSTAGGER3D: staggerflg out of bounds");
   if (istflg==1 && ox!=0) myerror("RPSTAGGER3D: ox should equal zero");
   if (istflg==2 && oy!=0) myerror("RPSTAGGER3D: oy should equal zero");
   if (istflg==3 && oz!=0) myerror("RPSTAGGER3D: oz should equal zero");

   if (staggerflg==1) {
      for (n=0;n<nx;n++) {
         np1=n+1;
         for (j=oy;j<=ny;j++) 
         for (k=oz;k<=nz;k++) 
            sf[n][j][k] = 0.5*(f[np1][j][k]+f[n][j][k]);
      }
      for (j=oy;j<=ny;j++) 
      for (k=oz;k<=nz;k++) sf[nx][j][k] = sf[0][j][k];

   } else if (staggerflg==2) {
      for (j=0;j<ny;j++) {
         jp1 = j+1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) 
            sf[n][j][k] = 0.5*(f[n][jp1][k]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (k=oz;k<=nz;k++) sf[n][ny][k] = sf[n][0][k];
      
   } else if (staggerflg==3) {
      for (k=0;k<nz;k++) {
         kp1=k+1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) 
            sf[n][j][k] = 0.5*(f[n][j][kp1]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (j=oy;j<=ny;j++) sf[n][j][nz] = sf[n][j][0];
      
   } else if (staggerflg==-1) {
      for (n=nx;n>=1;n--) {
         nm1=n-1;
         for (j=oy;j<=ny;j++) 
         for (k=oz;k<=nz;k++) 
            sf[n][j][k] = 0.5*(f[nm1][j][k]+f[n][j][k]);
      }
      for (j=oy;j<=ny;j++) 
      for (k=oz;k<=nz;k++) sf[0][j][k] = sf[nx][j][k];

   } else if (staggerflg==-2) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) 
            sf[n][j][k] = 0.5*(f[n][jm1][k]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (k=oz;k<=nz;k++) sf[n][0][k] = sf[n][ny][k];
      
   } else if (staggerflg==-3) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) 
            sf[n][j][k] = 0.5*(f[n][j][km1]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (j=oy;j<=ny;j++) sf[n][j][0] = sf[n][j][nz];
      
   }
}

/* ===========================================================================
   _RPSTAGGER3D: Jan 14/94

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  
   If staggerflg==1  average in x direction alone.
   If staggerflg==2  average in y direction alone.
   If staggerflg==3  average in z direction alone.

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required
   ===========================================================================
*/
void _RPstagger3D(f,ox,nx,oy,ny,oz,nz, staggerflg)
   int ox,nx,oy,ny,oz,nz,staggerflg;
   double ***f;
{
   int istflg,n,j,k;
   int np1,nm1,jp1,jm1,kp1,km1;
   double tmp;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>3) 
      myerror("_RPSTAGGER3D: staggerflg out of bounds");

   if (istflg==1 && ox!=0) myerror("_RPSTAGGER3D: ox should equal zero");
   if (istflg==2 && oy!=0) myerror("_RPSTAGGER3D: oy should equal zero");
   if (istflg==3 && oz!=0) myerror("_RPSTAGGER3D: oz should equal zero");

   if (staggerflg==1) {
      for (n=0;n<nx;n++) {
         np1=n+1;
         for (j=oy;j<=ny;j++) 
         for (k=oz;k<=nz;k++) 
            f[n][j][k] = 0.5*(f[np1][j][k]+f[n][j][k]);
      }
      for (j=oy;j<=ny;j++) 
      for (k=oz;k<=nz;k++) f[nx][j][k] = f[0][j][k];

   } else if (staggerflg==2) {
      for (j=0;j<ny;j++) {
         jp1 = j+1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) 
            f[n][j][k] = 0.5*(f[n][jp1][k]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (k=oz;k<=nz;k++) f[n][ny][k] = f[n][0][k];
      
   } else if (staggerflg==3) {
      for (k=0;k<nz;k++) {
         kp1=k+1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) 
            f[n][j][k] = 0.5*(f[n][j][kp1]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (j=oy;j<=ny;j++) f[n][j][nz] = f[n][j][0];
      
   } else if (staggerflg==-1) {
      for (n=nx;n>=1;n--) {
         nm1=n-1;
         for (j=oy;j<=ny;j++) 
         for (k=oz;k<=nz;k++) 
            f[n][j][k] = 0.5*(f[nm1][j][k]+f[n][j][k]);
      }
      for (j=oy;j<=ny;j++) 
      for (k=oz;k<=nz;k++) f[0][j][k] = f[nx][j][k];

   } else if (staggerflg==-2) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) 
            f[n][j][k] = 0.5*(f[n][jm1][k]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (k=oz;k<=nz;k++) f[n][0][k] = f[n][ny][k];
      
   } else if (staggerflg==-3) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) 
            f[n][j][k] = 0.5*(f[n][j][km1]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (j=oy;j<=ny;j++) f[n][j][0] = f[n][j][nz];
      
   }
}

/* ===========================================================================
   RZSTAGGER3D: Feb 21/94

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  
   If staggerflg==1  average in x and z directions.
   If staggerflg==2  average in y direction alone.
   If staggerflg==3  average in z direction alone.

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required

   Assume zero boundaries in direction of staggering.
   ===========================================================================
*/
void RZstagger3D(f,ox,nx,oy,ny,oz,nz, staggerflg,sf)
   int ox,nx,oy,ny,oz,nz,staggerflg;
   double ***f,***sf;
{
   int istflg,n,j,k;
   int np1,nm1,jp1,jm1,kp1,km1;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>3) 
      myerror("RZSTAGGER3D: staggerflg out of bounds");

   if (staggerflg==1 && ox!=1) myerror("RZSTAGGER3D: ox should equal one");
   if (staggerflg==2 && oy!=1) myerror("RZSTAGGER3D: oy should equal one");
   if (staggerflg==3 && oz!=1) myerror("RZSTAGGER3D: oz should equal one");
   if (staggerflg==-1 && ox!=0) myerror("RZSTAGGER3D: ox should equal zero");
   if (staggerflg==-2 && oy!=0) myerror("RZSTAGGER3D: oy should equal zero");
   if (staggerflg==-3 && oz!=0) myerror("RZSTAGGER3D: oz should equal zero");

   if (staggerflg==1) {
      for (n=1;n<nx;n++) {
         np1=n+1;
         for (j=oy;j<=ny;j++) 
         for (k=oz;k<=nz;k++) 
            sf[n][j][k] = 0.5*(f[np1][j][k]+f[n][j][k]);
      }
      for (j=oy;j<=ny;j++) 
      for (k=oz;k<=nz;k++) sf[0][j][k] = sf[nx][j][k] = 0.0;

   } else if (staggerflg==2) {
      for (j=1;j<ny;j++) {
         jp1 = j+1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) 
            sf[n][j][k] = 0.5*(f[n][jp1][k]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (k=oz;k<=nz;k++) sf[n][0][k] = sf[n][ny][k] = 0.0;

   } else if (staggerflg==3) {
      for (k=1;k<nz;k++) {
         kp1=k+1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) 
            sf[n][j][k] = 0.5*(f[n][j][kp1]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (j=oy;j<=ny;j++) sf[n][j][0] = sf[n][j][nz] = 0.0;
      
   } else if (staggerflg==-1) {
      for (n=nx;n>=1;n--) {
         nm1=n-1;
         for (j=oy;j<=ny;j++) 
         for (k=oz;k<=nz;k++) 
            sf[n][j][k] = 0.5*(f[nm1][j][k]+f[n][j][k]);
      }
   } else if (staggerflg==-2) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) 
            sf[n][j][k] = 0.5*(f[n][jm1][k]+f[n][j][k]);
      }
   } else if (staggerflg==-3) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) 
            sf[n][j][k] = 0.5*(f[n][j][km1]+f[n][j][k]);
      }
   }
}

/* ===========================================================================
   _RZSTAGGER3D: Feb 21/94

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  
   If staggerflg==1  average in x direction alone.
   If staggerflg==2  average in y direction alone.
   If staggerflg==3  average in z direction alone.

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required

   Assume zero boundaries in direction of staggering.
   ===========================================================================
*/
void _RZstagger3D(f,ox,nx,oy,ny,oz,nz, staggerflg)
   int ox,nx,oy,ny,oz,nz,staggerflg;
   double ***f;
{
   int istflg,n,j,k;
   int np1,nm1,jp1,jm1,kp1,km1;
   double tmp;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>3) 
      myerror("_RZSTAGGER3D: staggerflg out of bounds");

   if (staggerflg==1 && ox!=1) myerror("_RZSTAGGER3D: ox should equal one");
   if (staggerflg==2 && oy!=1) myerror("_RZSTAGGER3D: oy should equal one");
   if (staggerflg==3 && oz!=1) myerror("_RZSTAGGER3D: oz should equal one");
   if (staggerflg==-1 && ox!=0) myerror("_RZSTAGGER3D: ox should equal zero");
   if (staggerflg==-2 && oy!=0) myerror("_RZSTAGGER3D: oy should equal zero");
   if (staggerflg==-3 && oz!=0) myerror("_RZSTAGGER3D: oz should equal zero");

   if (staggerflg==1) {
      for (n=1;n<nx;n++) {
         np1=n+1;
         for (j=oy;j<=ny;j++) 
         for (k=oz;k<=nz;k++) 
            f[n][j][k] = 0.5*(f[np1][j][k]+f[n][j][k]);
      }
      for (j=oy;j<=ny;j++) 
      for (k=oz;k<=nz;k++) f[0][j][k] = f[nx][j][k] = 0.0;

   } else if (staggerflg==2) {
      for (j=1;j<ny;j++) {
         jp1 = j+1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) 
            f[n][j][k] = 0.5*(f[n][jp1][k]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (k=oz;k<=nz;k++) f[n][0][k] = f[n][ny][k] = 0.0;

   } else if (staggerflg==3) {
      for (k=1;k<nz;k++) {
         kp1=k+1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) 
            f[n][j][k] = 0.5*(f[n][j][kp1]+f[n][j][k]);
      }
      for (n=ox;n<=nx;n++)
      for (j=oy;j<=ny;j++) f[n][j][0] = f[n][j][nz] = 0.0;
      
   } else if (staggerflg==-1) {
      for (n=nx;n>=1;n--) {
         nm1=n-1;
         for (j=oy;j<=ny;j++) 
         for (k=oz;k<=nz;k++) 
            f[n][j][k] = 0.5*(f[nm1][j][k]+f[n][j][k]);
      }
   } else if (staggerflg==-2) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) 
            f[n][j][k] = 0.5*(f[n][jm1][k]+f[n][j][k]);
      }
   } else if (staggerflg==-3) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) 
            f[n][j][k] = 0.5*(f[n][j][km1]+f[n][j][k]);
      }
   }
}

/* ===========================================================================
   CSTAGGER3D: Oct 13/93

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  
   If staggerflg==1  average in x direction.
   If staggerflg==2  average in y direction.
   If staggerflg==3  average in z direction.

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required
   ===========================================================================
*/
void Cstagger3D(f,ox,nx,oy,ny,oz,nz, staggerflg,sf)
   int ox,nx,oy,ny,oz,nz,staggerflg;
   dcomplex ***f,***sf;
{
   int istflg,n,j,k;
   int np1,nm1,jp1,jm1,kp1,km1;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>3) 
      myerror("CSTAGGER3D: staggerflg out of bounds");

   if (staggerflg==1 && ox!=1) myerror("CSTAGGER3D: ox should equal one");
   if (staggerflg==2 && oy!=1) myerror("CSTAGGER3D: oy should equal one");
   if (staggerflg==3 && oz!=1) myerror("CSTAGGER3D: oz should equal one");
   if (staggerflg==-1 && ox!=0) myerror("CSTAGGER3D: ox should equal zero");
   if (staggerflg==-2 && oy!=0) myerror("CSTAGGER3D: oy should equal zero");
   if (staggerflg==-3 && oz!=0) myerror("CSTAGGER3D: oz should equal zero");

   if (staggerflg==1) {
      for (n=1;n<nx;n++) {
         np1=n+1;
         for (j=oy;j<=ny;j++)
         for (k=oz;k<=nz;k++) {
            sf[n][j][k].r = 0.5*(f[np1][j][k].r+f[n][j][k].r);
            sf[n][j][k].i = 0.5*(f[np1][j][k].i+f[n][j][k].i);
         }
      }
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) {
         sf[0][j][k].r = 1.5*f[1][j][k].r - 0.5*f[2][j][k].r;
         sf[0][j][k].i = 1.5*f[1][j][k].i - 0.5*f[2][j][k].i;
         sf[nx][j][k].r = 1.5*f[nx][j][k].r - 0.5*f[nx-1][j][k].r;
         sf[nx][j][k].i = 1.5*f[nx][j][k].i - 0.5*f[nx-1][j][k].i;
      }
   } else if (staggerflg==2) {
      for (j=1;j<ny;j++) {
         jp1 = j+1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) {
            sf[n][j][k].r = 0.5*(f[n][jp1][k].r+f[n][j][k].r);
            sf[n][j][k].i = 0.5*(f[n][jp1][k].i+f[n][j][k].i);
         }
      }
      for (n=ox;n<=nx;n++)
      for (k=oz;k<=nz;k++) {
         sf[n][0][k].r = 1.5*f[n][1][k].r - 0.5*f[n][2][k].r;
         sf[n][0][k].i = 1.5*f[n][1][k].i - 0.5*f[n][2][k].i;
         sf[n][ny][k].r = 1.5*f[n][ny][k].r - 0.5*f[n][ny-1][k].r;
         sf[n][ny][k].i = 1.5*f[n][ny][k].i - 0.5*f[n][ny-1][k].i;
      }
   } else if (staggerflg==3) {
      for (k=1;k<nz;k++) {
         kp1=k+1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) {
            sf[n][j][k].r = 0.5*(f[n][j][kp1].r+f[n][j][k].r);
            sf[n][j][k].i = 0.5*(f[n][j][kp1].i+f[n][j][k].i);
         }
      }
      for (n=ox;n<=nx;n++)
      for (j=oy;j<=ny;j++) {
         sf[n][j][0].r = 1.5*f[n][j][1].r - 0.5*f[n][j][2].r;
         sf[n][j][0].i = 1.5*f[n][j][1].i - 0.5*f[n][j][2].i;
         sf[n][j][nz].r = 1.5*f[n][j][nz].r - 0.5*f[n][j][nz-1].r;
         sf[n][j][nz].i = 1.5*f[n][j][nz].i - 0.5*f[n][j][nz-1].i;
      }
   } else if (staggerflg==-1) {
      for (n=nx;n>=1;n--) {
         nm1=n-1;
         for (j=oy;j<=ny;j++)
         for (k=oz;k<=nz;k++) {
            sf[n][j][k].r = 0.5*(f[nm1][j][k].r+f[n][j][k].r);
            sf[n][j][k].i = 0.5*(f[nm1][j][k].i+f[n][j][k].i);
         }
      }
   } else if (staggerflg==-2) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) {
            sf[n][j][k].r = 0.5*(f[n][jm1][k].r+f[n][j][k].r);
            sf[n][j][k].i = 0.5*(f[n][jm1][k].i+f[n][j][k].i);
         }
      }
   } else if (staggerflg==-3) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) {
            sf[n][j][k].r = 0.5*(f[n][j][km1].r+f[n][j][k].r);
            sf[n][j][k].i = 0.5*(f[n][j][km1].i+f[n][j][k].i);
         }
      }
   }
}

/* ===========================================================================
   _CSTAGGER3D: Sept 23/93

   Take a field defined at evenly spaced points j=0..ny, k=0..nz
   and replace with grid defined at intermediate points by average of
   surrounding points depending on value of staggerflg.  
   If staggerflg==1  average in x direction.
   If staggerflg==2  average in y direction.
   If staggerflg==3  average in z direction.

   Staggering from a (1,2,3,..n) grid to a (0,1,2,..,n) grid (ox==1)
   is done with staggerflg positive and extrapolation determines f[0],f[n]

   Staggering from a (0,1,2,..n) grid to a (1,2,3,..,n) grid (ox==0)
   is done with staggerflg negative and no extrapolation is required
   ===========================================================================
*/
void _Cstagger3D(f,ox,nx,oy,ny,oz,nz, staggerflg)
   int ox,nx,oy,ny,oz,nz,staggerflg;
   dcomplex ***f;
{
   int istflg,n,j,k;
   int np1,nm1,jp1,jm1,kp1,km1;
   dcomplex ctmp;

   istflg = iabs(staggerflg);

   if (istflg==0 || istflg>3) 
      myerror("_CSTAGGER3D: staggerflg out of bounds");

   if (staggerflg==1 && ox!=1) myerror("_CSTAGGER3D: ox should equal one");
   if (staggerflg==2 && oy!=1) myerror("_CSTAGGER3D: oy should equal one");
   if (staggerflg==3 && oz!=1) myerror("_CSTAGGER3D: oz should equal one");
   if (staggerflg==-1 && ox!=0) myerror("_CSTAGGER3D: ox should equal zero");
   if (staggerflg==-2 && oy!=0) myerror("_CSTAGGER3D: oy should equal zero");
   if (staggerflg==-3 && oz!=0) myerror("_CSTAGGER3D: oz should equal zero");

   if (staggerflg==1) {
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) {
         ctmp.r = 1.5*f[nx][j][k].r - 0.5*f[nx-1][j][k].r;
         ctmp.i = 1.5*f[nx][j][k].i - 0.5*f[nx-1][j][k].i;
         f[0][j][k].r = 1.5*f[1][j][k].r - 0.5*f[2][j][k].r;
         f[0][j][k].i = 1.5*f[1][j][k].i - 0.5*f[2][j][k].i;
         for (n=1;n<nx;n++) {
            f[n][j][k].r = 0.5*(f[n+1][j][k].r+f[n][j][k].r);
            f[n][j][k].i = 0.5*(f[n+1][j][k].i+f[n][j][k].i);
         }
         f[nx][j][k] = ctmp;
      }
   } else if (staggerflg==2) {
      for (n=ox;n<=nx;n++)
      for (k=oz;k<=nz;k++) {
         ctmp.r = 1.5*f[n][ny][k].r - 0.5*f[n][ny-1][k].r;
         ctmp.i = 1.5*f[n][ny][k].i - 0.5*f[n][ny-1][k].i;
         f[n][0][k].r = 1.5*f[n][1][k].r - 0.5*f[n][2][k].r;
         f[n][0][k].i = 1.5*f[n][1][k].i - 0.5*f[n][2][k].i;
         for (j=1;j<ny;j++) {
            f[n][j][k].r = 0.5*(f[n][j+1][k].r+f[n][j][k].r);
            f[n][j][k].i = 0.5*(f[n][j+1][k].i+f[n][j][k].i);
         }
         f[n][ny][k] = ctmp;
      }
   } else if (staggerflg==3) {
      for (n=ox;n<=nx;n++)
      for (j=oy;j<=ny;j++) {
         ctmp.r = 1.5*f[n][j][nz].r - 0.5*f[n][j][nz-1].r;
         ctmp.i = 1.5*f[n][j][nz].i - 0.5*f[n][j][nz-1].i;
         f[n][j][0].r = 1.5*f[n][j][1].r - 0.5*f[n][j][2].r;
         f[n][j][0].i = 1.5*f[n][j][1].i - 0.5*f[n][j][2].i;
         for (k=1;k<nz;k++) {
            f[n][j][k].r = 0.5*(f[n][j][k+1].r+f[n][j][k].r);
            f[n][j][k].i = 0.5*(f[n][j][k+1].i+f[n][j][k].i);
         }
      }
   } else if (staggerflg==-1) {
      for (n=nx;n>=1;n--) {
         nm1=n-1;
         for (j=oy;j<=ny;j++)
         for (k=oz;k<=nz;k++) {
            f[n][j][k].r = 0.5*(f[nm1][j][k].r+f[n][j][k].r);
            f[n][j][k].i = 0.5*(f[nm1][j][k].i+f[n][j][k].i);
         }
      }
   } else if (staggerflg==-2) {
      for (j=ny;j>=1;j--) {
         jm1 = j-1;
         for (n=ox;n<=nx;n++)
         for (k=oz;k<=nz;k++) {
            f[n][j][k].r = 0.5*(f[n][jm1][k].r+f[n][j][k].r);
            f[n][j][k].i = 0.5*(f[n][jm1][k].i+f[n][j][k].i);
         }
      }
   } else if (staggerflg==-3) {
      for (k=nz;k>=1;k--) {
         km1 = k-1;
         for (n=ox;n<=nx;n++)
         for (j=oy;j<=ny;j++) {
            f[n][j][k].r = 0.5*(f[n][j][km1].r+f[n][j][k].r);
            f[n][j][k].i = 0.5*(f[n][j][km1].i+f[n][j][k].i);
         }
      }
   }
}
