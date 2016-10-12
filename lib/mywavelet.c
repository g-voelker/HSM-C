#include <math.h>
#include "complex.h"
#include "alloc_space.h"
#include "constants.h"
#include "macros.h"

typedef struct {
   int ncof,ioff,joff;
   double *cc,*cr;
} wavefilt;

wavefilt wfilt;

void mypwtset(int n)
{
  int k;
  double sig = -1.0;
  static double c4[5]={0.0,0.4829629131445341,0.8365163037378079,
         0.2241438680420134,-0.1294095225512604};
  static double c12[13]={0.0,0.111540743350, 0.494623890398, 0.751133908021,
     0.315250351709,-0.226264693965,-0.129766867567,
     0.097501605587, 0.027522865530,-0.031582039318,
     0.000553842201, 0.004777257511,-0.001077301085};
  static double c20[21]={0.0,0.026670057901, 0.188176800078, 0.527201188932,
     0.688459039454, 0.281172343661,-0.249846424327,
     -0.195946274377, 0.127369340336, 0.093057364604,
     -0.071394147166,-0.029457536822, 0.033212674059,
     0.003606553567,-0.010733175483, 0.001395351747,
     0.001992405295,-0.000685856695,-0.000116466855,
     0.000093588670,-0.000013264203};
  static double c4r[5],c12r[13],c20r[21];

  wfilt.ncof=n;
  if (n==4) {
     wfilt.cc=c4;
     wfilt.cr=c4r;
  }
  if (n==12) {
     wfilt.cc=c12;
     wfilt.cr=c12r;
  }
  if (n==20) {
     wfilt.cc=c20;
     wfilt.cr=c20r;
  }
  else nrerror("unimplemented value n in pwtset");
  for (k=1;k<=n;k++) {
     wfilt.cr[wfilt.ncof+1-k]=sig*wfilt.cc[k];
     sig = -sig;
  }
  wfilt.ioff = wfilt.joff = -(n >> 1);
}


void mypwt(double a[], unsigned long n, int isign)
{
  int it;
  double ai,ai1,*wksp;
  unsigned long i,ii,j,jf,jr,k,n1,ni,nj,nh,nmod;

  if (n<4) return;
  wksp=dvector(1,n);
  nmod=wfilt.ncof*n;
  n1=n-1;
  nh=n >> 1;
  for (j=1; j<=n; j++) wksp[j]=0.0;
  if (isign >= 0) {
     for (ii=1,i=1;i<=n;i+=2,ii++) {
        ni=i+nmod+wfilt.ioff;
        nj=i+nmod+wfilt.joff;
        for (k=1;k<=wfilt.ncof;k++) {
           jf=n1 & (ni+k);
           jr=n1 & (nj+k);
           wksp[ii] += wfilt.cc[k]*a[jf+1];
           wksp[ii+nh] += wfilt.cr[k]*a[jr+1];
        }
     }
  } else {
     for (ii=1,i=1;i<=n;i+=2,ii++) {
        ai=a[ii];
        ai1=a[ii+nh];
        ni=i+nmod+wfilt.ioff;
        nj=i+nmod+wfilt.joff;
        for (k=1;k<=wfilt.ncof;k++) {
           jf=(n1 & (ni+k))+1;
           jr=(n1 & (nj+k))+1;
           wksp[jf] += wfilt.cc[k]*ai;
           wksp[jr] += wfilt.cr[k]*ai1;
        }
     }
  }
  for (j=1;j<=n;j++) a[j]=wksp[j];
  free_dvector(wksp,1,n);
}

/* ========================================================================
   WAVELET: Mar 2/97

   Wavelet transform of a vector (1..nn) of double data (see routine
   "wt1" in Numerical recipes).  Assumed Daubechies wavelet filter with 
   20 coefficients.  nn is a power of 2.
   ========================================================================
*/
void wavelet(double *a, unsigned long nn, double *wla, int isign)
{
   unsigned long n;
 
   if (nn<4) return;

   if (isign >=0) {  /* a->wla */
     for (n=1;n<=nn;n++) wla[n]=a[n];
   } else {          /* wla->a */
     for (n=1;n<=nn;n++) a[n]=wla[n];
   }

   /* initialize wavelet coefficients (wfilt) */
   mypwtset(20);

   if (isign >=0) {  /* transform */
     for (n=nn;n>=4;n>>=1) mypwt(wla,n,isign);
   } else { 
     for (n=4;n<=nn;n<<=1) mypwt(a,n,isign);
   }
}


/* ========================================================================
   WAVELET2D: Mar 2/97

   Wavelet filter of a matrix vector (1..nn,1..mm) of double data
   ========================================================================
*/
void wavelet2D(
   double **mat, unsigned long nn, unsigned long mm, double **wlmat, int isign)
{
   unsigned long n,m;
   double *a;

   a = dvector(1,nn);

   if (isign==1) {
     for (n=1;n<=nn;n++) for (m=1;m<=mm;m++) wlmat[n][m]=mat[n][m];
   } else {
     for (n=1;n<=nn;n++) for (m=1;m<=mm;m++) mat[n][m]=wlmat[n][m];
   }

   /* initialize wavelet coefficients (wfilt) */
   mypwtset(20);

   for (n=1;n<=nn;n++) { 
     if (mm<4) return;
     if (isign >=0) {  /* transform */
       for (m=mm;m>=4;m>>=1) mypwt(wlmat[n],m,isign);
     } else { 
       for (m=4;m<=mm;m<<=1) mypwt(mat[n],m,isign);
     }
   }

   for (m=1;m<=mm;m++) { 
     if (nn<4) return;
     if (isign >=0) {  /* transform */
       for (n=1;n<=nn;n++) a[n]=wlmat[n][m];
       for (n=nn;n>=4;n>>=1) mypwt(a,n,isign);
       for (n=1;n<=nn;n++) wlmat[n][m]=a[n];
     } else { 
       for (n=1;n<=nn;n++) a[n]=mat[n][m];
       for (n=4;n<=nn;n<<=1) mypwt(a,n,isign);
       for (n=1;n<=nn;n++) mat[n][m]=a[n];
     }
   }

   free_dvector(a,1,nn);
}
