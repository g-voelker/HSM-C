/* ===================================================================
   Routines to take fast fourier transforms of single vectors of data. 
   These routines are similar to Numerical Recipes except that the
   backward tranform is defined by:

       f_k = (1/N) sum( n=0..N-1, F_n exp(I n (alpha x)) ) 

   and the forward transform is defined by:

       F_n = sum( k=0..N-1, f_k exp(-I n (alpha x_k)) ) 

   where  alpha=2*PI/L, x_k = L*k/N  so   alpha x_k = 2*PI*k/N.
   (In Numerical Recipes the sign before the complex I is reversed.
    The normalization is unchanged but should be set in utilizing
    programs so that 1/N precedes the forward transform).

   UPDATE AUG 6/93:  The 2ndEdition does things slightly differently
   so I've switch the signs back to the values in the book and made a
   a few other changes so this program and the new ones are more alike
   ===================================================================
 */
#include <math.h>
#include "complex.h"
#include "alloc_space.h"
#include "constants.h"
#include "macros.h"

static double rarg1;
#define SWAP(a,b) rarg1=(a); (a)=(b); (b)=rarg1

void myfour1(data,nn,isign)
   int nn,isign;
   double *data;
{
   int n, mmax, m, j, istep, i;
   double wtemp, wr, wpr, wpi, wi, theta;
   double tempr, tempi;

   n = nn << 1;
   j = 1;

/* bit reversal */

   for (i=1; i<n; i+=2) {
      if (j>i) {
         SWAP( data[j], data[i] );
         SWAP( data[j+1], data[i+1] );
      }
      m=n >> 1;
      while (m >= 2 && j>m) {
         j -= m;
         m >>= 1;
      }
      j += m;
   }

/* fourier breakdown */

   mmax = 2;
   while (n > mmax) {
      istep = 2*mmax;
      theta = isign*(PI2/mmax); 
      wtemp = sin( 0.5*theta );
      wpr = -2.0*wtemp*wtemp;
      wpi = sin( theta );
      wr = 1.0;
      wi = 0.0;
      for (m=1; m<mmax; m +=2) {
         for (i=m; i<=n; i += istep) {
            j = i+mmax;
            tempr = wr*data[j] - wi*data[j+1];
            tempi = wr*data[j+1] + wi*data[j];
            data[j] = data[i] - tempr;
            data[j+1] = data[i+1] - tempi;
            data[i] += tempr;
            data[i+1] += tempi;
         }
         wr = (wtemp=wr)*wpr - wi*wpi + wr;
         wi = wi*wpr + wtemp*wpi + wi;
      }
      mmax = istep;
   }
}

void myrealft(data,n,isign)
   int n,isign;
   double *data;
{
   int i, i1, i2, i3, i4, np3;
   double wtemp, wr, wpr, wpi, wi, theta;
   double c1=0.5, c2, h1r, h1i, h2r, h2i;

   theta = PI/(double) (n>>1);
   if (isign==1) {        /* do forward transform */
      c2 = -0.5;
      myfour1(data,n>>1,1);
   } else {               /* set up for backward transform */
      c2 = 0.5;
      theta = -theta;
   }

   wtemp = sin(0.5*theta);
   wpr = -2.0*wtemp*wtemp;
   wpi = sin(theta);
   wr = 1.0 + wpr;
   wi = wpi;
   np3 = n + 3;
   for (i=2; i<=(n>>2); i++) {
      i4 = 1 + (i3 = np3-(i2 = 1+(i1 = i+i-1)));
      h1r = c1*( data[i1] + data[i3] );
      h1i = c1*( data[i2] - data[i4] );
      h2r = -c2*( data[i2] + data[i4] );
      h2i = c2*( data[i1] - data[i3] );
      data[i1] = h1r + wr*h2r - wi*h2i;
      data[i2] = h1i + wr*h2i + wi*h2r;
      data[i3] = h1r - wr*h2r + wi*h2i;
      data[i4] = - h1i + wr*h2i + wi*h2r;
      wr = (wtemp=wr)*wpr - wi*wpi + wr;
      wi = wi*wpr + wtemp*wpi +wi;
   }
   if (isign == 1) {
      data[1] = (h1r=data[1])+data[2];
      data[2] = h1r - data[2];
   } else {
      data[1] = c1*((h1r=data[1])+data[2]);
      data[2] = c1*(h1r-data[2]);
      myfour1(data,n>>1,-1);
   }
}

/* ====================================================================
   MYCOSFT 
   ====================================================================
*/
void mycosft(y,n)
   int n;
   double *y;
{
   int j,n2;
   double sum,y1,y2;
   double theta,wi=0.0,wr=1.0,wpi,wpr,wtemp;

   theta= PI/(double) n; 
   wtemp=sin(0.5*theta);
   wpr = -2.0*wtemp*wtemp;
   wpi=sin(theta);
   sum=0.5*(y[1]-y[n+1]);
   y[1]=0.5*(y[1]+y[n+1]);
   n2=n+2;
   for (j=2;j<=(n>>1);j++) {
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
      y1=0.5*(y[j]+y[n2-j]);
      y2=(y[j]-y[n2-j]);
      y[j]=y1-wi*y2;
      y[n2-j]=y1+wi*y2;
      sum += wr*y2;
   }
   myrealft(y,n,1);
   y[n+1]=y[2];
   y[2]=sum;
   for (j=4;j<=n;j+=2) {
      sum += y[j];
      y[j]=sum;
   }
}

/* ====================================================================
   MYCOSFT2 is good for data defined midway between meshpoints 
   ====================================================================
*/
void mycosft2(double *y, int n, int isign)
{
   int i;
   double sum,sum1,y1,y2,ytemp;
   double theta,wi=0.0,wi1,wpi,wpr,wr=1.0,wr1,wtemp;

   theta=0.5*PI/n;
   wr1=cos(theta);
   wi1=sin(theta);
   wpr = -2.0*wi1*wi1;
   wpi=sin(2.0*theta);
   if (isign == 1) {
      for (i=1;i<=n/2;i++) {
         y1=0.5*(y[i]+y[n-i+1]);
         y2=wi1*(y[i]-y[n-i+1]);
         y[i]=y1+y2;
         y[n-i+1]=y1-y2;
         wr1=(wtemp=wr1)*wpr-wi1*wpi+wr1;
         wi1=wi1*wpr+wtemp*wpi+wi1;
      }
      myrealft(y,n,1);
      for (i=3;i<=n;i+=2) {
         wr=(wtemp=wr)*wpr-wi*wpi+wr;
         wi=wi*wpr+wtemp*wpi+wi;
         y1=y[i]*wr-y[i+1]*wi;
         y2=y[i+1]*wr+y[i]*wi;
         y[i]=y1;
         y[i+1]=y2;
      }
      sum=0.5*y[2];
      for (i=n;i>=2;i-=2) {
         sum1=sum;
         sum += y[i];
         y[i]=sum1;
      }
   } else if (isign == -1) {
      ytemp=y[n];
      for (i=n;i>=4;i-=2) y[i]=y[i-2]-y[i];
      y[2]=2.0*ytemp;
      for (i=3;i<=n;i+=2) {
         wr=(wtemp=wr)*wpr-wi*wpi+wr;
         wi=wi*wpr+wtemp*wpi+wi;
         y1=y[i]*wr+y[i+1]*wi;
         y2=y[i+1]*wr-y[i]*wi;
         y[i]=y1;
         y[i+1]=y2;
      }
      myrealft(y,n,-1);
      for (i=1;i<=n/2;i++) {
         y1=y[i]+y[n-i+1];
         y2=(0.5/wi1)*(y[i]-y[n-i+1]);
         y[i]=0.5*(y1+y2);
         y[n-i+1]=0.5*(y1-y2);
         wr1=(wtemp=wr1)*wpr-wi1*wpi+wr1;
         wi1=wi1*wpr+wtemp*wpi+wi1;
      }
   }
}

/* ====================================================================
   MYSINFT
   ====================================================================
*/
void mysinft(y,n)
   int n;
   double *y;
{
   int j,n2=n+2;
   double sum,y1,y2;
   double theta,wi=0.0,wr=1.0,wpi,wpr,wtemp;

   theta= PI/(double) n;   
   wtemp=sin(0.5*theta);
   wpr = -2.0*wtemp*wtemp;
   wpi=sin(theta);
   y[1]=0.0;
   for (j=2;j<=(n>>1)+1;j++) {
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
      y1=wi*(y[j]+y[n2-j]);
      y2=0.5*(y[j]-y[n2-j]);
      y[j]=y1+y2;
      y[n2-j]=y1-y2;
   }
   myrealft(y,n,1);
   y[1]*=0.5;
   sum=y[2]=0.0;
   for (j=1;j<=n-1;j+=2) {
      sum += y[j];
      y[j]=y[j+1];
      y[j+1]=sum;
   }
}

void mytwofft(data1,data2,fft1,fft2, n)
   int n;
   double *data1,*data2,*fft1,*fft2;
{
   int nn3,nn2,jj,j;
   double rep,rem,aip,aim;

   nn3=1+(nn2=2+n+n);
   for (j=1,jj=2;j<=n;j++,jj+=2) {
      fft1[jj-1]=data1[j];
      fft1[jj]=data2[j];
   }
   myfour1(fft1,n,1);
   fft2[1]=fft1[2];
   fft1[2]=fft2[2]=0.0;
   for (j=3;j<=n+1;j+=2) {
      rep=0.5*(fft1[j]+fft1[nn2-j]);
      rem=0.5*(fft1[j]-fft1[nn2-j]);
      aip=0.5*(fft1[j+1]+fft1[nn3-j]);
      aim=0.5*(fft1[j+1]-fft1[nn3-j]);
      fft1[j]=rep;
      fft1[j+1]=aim;
      fft1[nn2-j]=rep;
      fft1[nn3-j] = -aim;
      fft2[j]=aip;
      fft2[j+1] = -rem;
      fft2[nn2-j]=aip;
      fft2[nn3-j]=rem;
   }
}

void myconvlv(data,n,respns,m,rfflg,isign, ans)
   int n,m,rfflg,isign;
   double *data,*respns,*ans;
{
   int i,no2;
   double dum,mag2,*fft;

   fft=dvector(1,n<<1);
   for (i=1;i<=(m-1)/2;i++)
      respns[n+1-i]=respns[m+1-i];
   for (i=(m+3)/2;i<=n-(m-1)/2;i++)
      respns[i]=0.0;
   mytwofft(data,respns,fft,ans,n);
   no2=n>>1;
   for (i=2;i<=n+2;i+=2) {
      if (isign == 1) {
         ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
         ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
      } else if (isign == -1) {
               if ((mag2=sqr(ans[i-1])+sqr(ans[i])) == 0.0)
            myerror("MYCONVLV: Deconvolving at response zero in CONVLV");
         ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
         ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
      } else myerror("MYCONVLV: No meaning for ISIGN in CONVLV");
   }
   ans[2]=ans[n+1];
   myrealft(ans,no2,-1);
   free_dvector(fft,1,n<<1);
}

void mycorrel(data1,data2,n, ans)
   int n;
   double *data1,*data2,*ans;
{
   int no2,i;
   double dum,*fft;

   fft=dvector(1,n<<1);
   mytwofft(data1,data2,fft,ans,n);
   no2=n>>1;
   for (i=2;i<=n+2;i+=2) {
      ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2;
      ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2;
   }
   ans[2]=ans[n+1];
   myrealft(ans,no2,-1);
   free_dvector(fft,1,n<<1);
}
#undef SWAP

/* ======================================================================
   SUBROUTINE: myfft
   Calculate fourier transform of real vector to spectral vector (sgn=1)
      and vice versa (sgn=-1).  Note - for a real field the fourier 
      transform has negative frequency components which are complex 
      conjugates of the positive frequency parts.  Only the positive 
      frequency components are returned.

   vec is a nn real vector ( vec[0]==vec[nn] )
   cvec is an (nn/2) complex vector  (cvec[0].i=cvec[nn/2].i=0)

   The point vec[0] corresponds to x=0.
   The element cvec[0] corresponds to the zero frequency mode.

   According to Numerical Recipes, fourier transforms of 2*nhf data to
      nhf complex frequency components are not normalized but inverse
      transform is normalized by  nhf.  In this case the zero frequency
      term will correspond to the sum of all real data (NOT THE AVERAGE!).
      Well, this sucks.  It makes much more sense to take the average by
      dividing the forward transform by nhf.  So that's what I do.  Nyah.
   ======================================================================
*/
void myfft(double *vec, int nn, dcomplex *cvec, int isgn)
{
   int n,n2,n2p1,nhf;
   double norm;
   double *ft;

   nhf = nn>>1;
   norm = 1.0/((double) nn);      /* see comments above */

   /* allocate space */
   ft = dvector(0,nn);

   /* Complex FT values are stored in vectors with components 2n-1 (real)
      and 2n (imaginary).  In frequency space, a[1] is zero frequency
      component (the sum of the real space values) and a[2] is the
      1/2d frequency component (alternating sum of real space values).
   */
   if (isgn==1) {
      for (n=1;n<=nn;n++) ft[nn-n]=vec[n];

      myrealft(ft-1,nn,1);

      cvec[0].r = norm*ft[0];
      cvec[nhf].r = norm*ft[1];
      cvec[0].i = 0.0;
      cvec[nhf].i = 0.0;
      for (n=1;n<nhf;n++) {
         n2 = n<<1;
         n2p1 = n2+1;
         cvec[n].r = norm*ft[n2];
         cvec[n].i = norm*ft[n2p1];
      }
   } else {
      ft[0] = cvec[0].r;
      ft[1] = cvec[nhf].r;
      for (n=1;n<nhf;n++) {
         n2 = n<<1;
         n2p1 = n2+1;
         ft[n2] = cvec[n].r;
         ft[n2p1] = cvec[n].i;
      }

      /* fft ftvec to physical space */
      myrealft(ft-1,nn,-1);
      for (n=1;n<=nn;n++) vec[n] = 2*ft[nn-n];

      /* reproduce first column in last column to show periodicity */
      vec[0]=vec[nn];
   }

   free_dvector(ft,0,nn);
}

/* ======================================================================
   SUBROUTINE: myzfft
   zomplex version of myfft
   ======================================================================
*/
void myzfft(double *vec, int nn, zomplex *cvec, int isgn)
{
   int n,n2,n2p1,nhf;
   double norm;
   double *ft;

   nhf = nn>>1;
   norm = 1.0/((double) nn);      /* see comments above */

   /* allocate space */
   ft = dvector(0,nn);

   /* Complex FT values are stored in vectors with components 2n-1 (real)
      and 2n (imaginary).  In frequency space, a[1] is zero frequency
      component (the sum of the real space values) and a[2] is the
      1/2d frequency component (alternating sum of real space values).
   */
   if (isgn==1) {
      for (n=1;n<=nn;n++) ft[nn-n]=vec[n];

      myrealft(ft-1,nn,1);

      cvec[0].re = norm*ft[0];
      cvec[nhf].re = norm*ft[1];
      cvec[0].im = 0.0;
      cvec[nhf].im = 0.0;
      for (n=1;n<nhf;n++) {
         n2 = n<<1;
         n2p1 = n2+1;
         cvec[n].re = norm*ft[n2];
         cvec[n].im = norm*ft[n2p1];
      }
   } else {
      ft[0] = cvec[0].re;
      ft[1] = cvec[nhf].re;
      for (n=1;n<nhf;n++) {
         n2 = n<<1;
         n2p1 = n2+1;
         ft[n2] = cvec[n].re;
         ft[n2p1] = cvec[n].im;
      }

      /* fft ftvec to physical space */
      myrealft(ft-1,nn,-1);
      for (n=1;n<=nn;n++) vec[n] = 2*ft[nn-n];

      /* reproduce first column in last column to show periodicity */
      vec[0]=vec[nn];
   }

   free_dvector(ft,0,nn);
}

/* ======================================================================
   SUBROUTINE: mycfft
   Calculate fourier transform of complex vector to spectral vector (sgn=1)
      and vice versa (sgn=-1).  

   cv is an 0..mm complex vector  (0 value  = mm value)
   fcv is an 0..mm complex vector (mm value = 0 value)
   mhf corresponds to x=0.
   ======================================================================
*/
void mycfft(dcomplex *cv, int mm, dcomplex *fcv, int isgn)
{
   int m,m2,mhf,mm2,mpmhf,mmmhf,mmm;
   double norm;
   double *rfld;

   mm2 = mm<<1;
   mhf = mm>>1;
   norm = 1.0/((double) mm);

   rfld=dvector(0,mm2);

   if (isgn==1) {
     for (m=0;m<mm;m++) {
       m2 = m<<1;
       mmm=mm-m;
       rfld[m2] = cv[mmm].r;
       rfld[m2+1] = cv[mmm].i;
     }

     myfour1(rfld-1,mm,1);

     for (m=0;m<mhf;m++) {
       m2=m<<1;
       mpmhf=mhf+m;
       fcv[mpmhf].r = norm*rfld[m2];
       fcv[mpmhf].i = norm*rfld[m2+1];
     }
     for (m=mhf;m<mm;m++) {
       m2=m<<1;
       mmmhf= m-mhf;
       fcv[mmmhf].r = norm*rfld[m2];
       fcv[mmmhf].i = norm*rfld[m2+1];
     }
     fcv[mm]=fcv[0];

   } else if (isgn==-1) {

     for (m=0;m<mhf;m++) {
        m2=m<<1;
        mpmhf=m+mhf;
        rfld[m2] = fcv[mpmhf].r;
        rfld[m2+1] = fcv[mpmhf].i;
     }
     for (m=mhf;m<mm;m++) {
        m2=m<<1;
        mmmhf=m-mhf;
        rfld[m2] = fcv[mmmhf].r;
        rfld[m2+1] = fcv[mmmhf].i;
     }

     myfour1(rfld-1,mm,-1);

     for (m=0;m<mm;m++) {
        m2=m<<1;
        mmm=mm-m;
        cv[mmm].r = rfld[m2];
        cv[mmm].i = rfld[m2+1];
     }
     cv[0] = cv[mm];
   }

   free_dvector(rfld,0,mm2);
}

/* ======================================================================
   SUBROUTINE: myczfft
   zomplex version of mycfft 
   ======================================================================
*/
void myczfft(zomplex *cv, int mm, zomplex *fcv, int isgn)
{
   int m,m2,mhf,mm2,mpmhf,mmmhf,mmm;
   double norm;
   double *rfld;

   mm2 = mm<<1;
   mhf = mm>>1;
   norm = 1.0/((double) mm);

   rfld=dvector(0,mm2);

   if (isgn==1) {
     for (m=0;m<mm;m++) {
       m2 = m<<1;
       mmm=mm-m;
       rfld[m2] = cv[mmm].re;
       rfld[m2+1] = cv[mmm].im;
     }

     myfour1(rfld-1,mm,1);

     for (m=0;m<mhf;m++) {
       m2=m<<1;
       mpmhf=mhf+m;
       fcv[mpmhf].re = norm*rfld[m2];
       fcv[mpmhf].im = norm*rfld[m2+1];
     }
     for (m=mhf;m<mm;m++) {
       m2=m<<1;
       mmmhf= m-mhf;
       fcv[mmmhf].re = norm*rfld[m2];
       fcv[mmmhf].im = norm*rfld[m2+1];
     }
     fcv[mm]=fcv[0];

   } else if (isgn==-1) {

     for (m=0;m<mhf;m++) {
        m2=m<<1;
        mpmhf=m+mhf;
        rfld[m2] = fcv[mpmhf].re;
        rfld[m2+1] = fcv[mpmhf].im;
     }
     for (m=mhf;m<mm;m++) {
        m2=m<<1;
        mmmhf=m-mhf;
        rfld[m2] = fcv[mmmhf].re;
        rfld[m2+1] = fcv[mmmhf].im;
     }

     myfour1(rfld-1,mm,-1);

     for (m=0;m<mm;m++) {
        m2=m<<1;
        mmm=mm-m;
        cv[mmm].re = rfld[m2];
        cv[mmm].im = rfld[m2+1];
     }
     cv[0] = cv[mm];
   }

   free_dvector(rfld,0,mm2);
}
