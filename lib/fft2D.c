#include <math.h>
#include "complex.h"
#include "alloc_space.h"
#include "myerror.h"
#include "myfft.h"

/* ======================================================================
   SUBROUTINE: hfft2D
   Calculate fourier transform of real field to horizontally spectral,
      vertically real field.  If sgn=-1 find inverse FT from horizontally
      spectral, vertically real field to real field.  Note: for a real
      field the fourier transform has negative frequency components 
      which are complex conjugates of the positive frequency parts.  Only
      the positive frequency components are returned.

   mat is a nn by mm real matrix ( mat[0][m]==mat[nn][m] )
   cfld is an (nn/2) by mm complex matrix  (cfld[0][mm].i=cfld[nn/2][mm].i=0)
  
   The point mat[0][0] corresponds to x=z=0.
   The element cfld[0][m] corresponds to the zero frequency mode at level m.

   SEE ALSO ZNFFT2D BELOW

   According to Numerical Recipes, fourier transforms of 2*nhf data to
      nhf complex frequency components are not normalized but inverse
      transform is normalized by  nhf.  In this case the zero frequency
      term will correspond to the sum of all real data (NOT THE AVERAGE!).
      Well, this sucks.  It makes much more sense to take the average by
      dividing the forward transform by nhf.  So that's what I do.  Nyah. 
   ======================================================================
*/
void hfft2D(mat,nn,mm, cfld,isgn)
   int nn,mm,isgn;
   double **mat;
   dcomplex **cfld;
{
   int n,m,n2,n2p1,nhf;
   double norm;
   double *ft;

   nhf = nn>>1;
   norm = 1.0/((double) nn);      /* see comments above */

   /* allocate space */
   ft = dvector(0,nn);

   /* Complex FT values are stored in matrices with components 2n-1 (real)
      and 2n (imaginary).  In frequency space, a[1] is zero frequency
      component (the sum of the real space values) and a[2] is the
      1/2d frequency component (alternating sum of real space values).
   */
   if (isgn==1) {
      for (m=0;m<mm;m++) {    /* cfld at mm equals cfld at zero */
         for (n=1;n<=nn;n++) ft[nn-n]=mat[n][m];

         myrealft(ft-1,nn,1);

         cfld[0][m].r = norm*ft[0];
         cfld[nhf][m].r = norm*ft[1];
         cfld[0][m].i = 0.0;
         cfld[nhf][m].i = 0.0;
         for (n=1;n<nhf;n++) {
            n2 = n<<1;
            n2p1 = n2+1;
            cfld[n][m].r = norm*ft[n2];
            cfld[n][m].i = norm*ft[n2p1];
         }
      }
      for (n=0;n<=nhf;n++) cfld[n][mm]=cfld[n][0];
   } else {
      for (m=0;m<=mm;m++) {  
         ft[0] = cfld[0][m].r;
         ft[1] = cfld[nhf][m].r;
         for (n=1;n<nhf;n++) {
            n2 = n<<1;
            n2p1 = n2+1;
            ft[n2] = cfld[n][m].r;
            ft[n2p1] = cfld[n][m].i;
         }

         /* fft ftvec to physical space */
         myrealft(ft-1,nn,-1);
         for (n=1;n<=nn;n++) mat[n][m] = 2*ft[nn-n];

         /* reproduce first column in last column to show periodicity */
         mat[0][m]=mat[nn][m];
      }
   }

   free_dvector(ft,0,nn);
}

/* =====================================================================
   SUBROUTINE: VRFFT2D
   Calculate vertical fourier transform of real field.  
      If sgn=-1 find inverse FT from vertically spectral to
      vertically real field

   mat is a nn by mm real matrix ( mat[n][0]==mat[n][mm] )
   cfld is an nn by (mm/2) complex matrix  (cfld[nn][0].i=cfld[nn][mm/2].i=0)
  
   The point mat[0][0] corresponds to x=z=0.
   The element cfld[n][0] corresponds to the zero frequency mode at level n.
   =====================================================================
*/
void vrfft2D(double **mat,int nn,int mm,dcomplex **cfld, int isgn)
{
   int n,m,m2,m2p1,mhf;
   double norm;
   double *ft;

   mhf = mm>>1;
   norm = 1.0/((double) mm);      /* see to HFFT2D comments */

   /* allocate space */
   ft = dvector(0,mm);

   /* Complex FT values are stored in matrices with components 2n-1 (real)
      and 2n (imaginary).  In frequency space, a[1] is zero frequency
      component (the sum of the real space values) and a[2] is the
      1/2d frequency component (alternating sum of real space values).
   */
   if (isgn==1) {
      for (n=0;n<nn;n++) {  /* cfld at nn equals cfld at zero */
         for (m=1;m<=mm;m++) ft[mm-m]=mat[n][m];

         myrealft(ft-1,mm,1);

         cfld[n][0].r = norm*ft[0];
         cfld[n][mhf].r = norm*ft[1];
         cfld[n][0].i = 0.0;
         cfld[n][mhf].i = 0.0;
         for (m=1;m<mhf;m++) {
            m2 = m<<1;
            m2p1 = m2+1;
            cfld[n][m].r = norm*ft[m2];
            cfld[n][m].i = norm*ft[m2p1];
         }
      }
      for (m=0;m<=mhf;m++) cfld[nn][m]=cfld[0][m];
   } else {
      for (n=0;n<=nn;n++) { 
         ft[0] = cfld[n][0].r;
         ft[1] = cfld[n][mhf].r;
         for (m=1;m<mhf;m++) {
            m2 = m<<1;
            m2p1 = m2+1;
            ft[m2] = cfld[n][m].r;
            ft[m2p1] = cfld[n][m].i;
         }

         /* fft ftvec to physical space */
         myrealft(ft-1,mm,-1);
         for (m=1;m<=mm;m++) mat[n][m] = 2*ft[mm-m];

         /* reproduce first column in last column to show periodicity */
         mat[n][0]=mat[n][mm];
      }
   }
   free_dvector(ft,0,mm);
}

/* =====================================================================
   SUBROUTINE: VFFT2D
   Calculate fourier transform of horizontally spectral, vertically real
      fields.  If sgn=-1 find inverse FT from 2D-spectra to horizontally
      spectral, vertically real fields

   forward transform from fld->ffld
   backward transform from ffld->fld

   The element fld[0][m] corresponds to the zero horiz. freq at level m.
   The element ffld[0][mhf] corresponds to the zero frequency mode

   both fld,ffld have size 0..nhf, 0..mm
   =====================================================================
*/
void vfft2D(fld,nhf,mm, ffld,isgn)
   int nhf,mm,isgn;
   dcomplex **fld,**ffld;
{
   int mhf,n,m,m2,mpmhf,mmmhf,mmm;
   double norm;
   double *rfld;

   rfld=dvector(0,2*mm);

   mhf = mm>>1;
   norm = 1.0/((double) mm);     /* see comments in header of hfft2D */

   if (isgn==1) {
      for (n=0;n<=nhf;n++) {
         for (m=0;m<mm;m++) {
            m2 = m<<1;
            mmm=mm-m;
            rfld[m2] = fld[n][mmm].r;
            rfld[m2+1] = fld[n][mmm].i;
         }

         myfour1(rfld-1,mm,1);

         for (m=0;m<mhf;m++) {
            m2=m<<1;
            mpmhf=mhf+m;
            ffld[n][mpmhf].r = norm*rfld[m2];
            ffld[n][mpmhf].i = norm*rfld[m2+1];
         }
         for (m=mhf;m<mm;m++) {
            m2=m<<1;
            mmmhf= m-mhf;
            ffld[n][mmmhf].r = norm*rfld[m2];
            ffld[n][mmmhf].i = norm*rfld[m2+1];
         }
         /* enforce zero values of im parts for real vectors */
         if (n==0 || n==nhf) ffld[n][mhf].i = ffld[n][0].i = 0.0;

         ffld[n][mm]=ffld[n][0];
      }
   } else if (isgn==-1) {
      for (n=0;n<=nhf;n++) {
         for (m=0;m<mhf;m++) {
            m2=m<<1;
            mpmhf=m+mhf;
            rfld[m2] = ffld[n][mpmhf].r;
            rfld[m2+1] = ffld[n][mpmhf].i;
         }
         for (m=mhf;m<mm;m++) {
            m2=m<<1;
            mmmhf=m-mhf;
            rfld[m2] = ffld[n][mmmhf].r;
            rfld[m2+1] = ffld[n][mmmhf].i;
         }

         myfour1(rfld-1,mm,-1);

         for (m=0;m<mm;m++) {
            m2=m<<1;
            mmm=mm-m;
            fld[n][mmm].r = rfld[m2];
            fld[n][mmm].i = rfld[m2+1];
         }
         /* enforce zero values of im parts for real vectors */
         if (n==0 || n==nhf) for(m=0;m<=mm;m++) fld[n][m].i = 0.0;

         fld[n][0] = fld[n][mm];
      }
   }

   free_dvector(rfld,0,2*mm);
}

/* =====================================================================

   ****  WARNING WARNING WARNING  - THIS HAS NOT BEEN TESTED !!!!!!!

   SUBROUTINE: VSFFT2D
   Calculate fourier SINE transform of horizontally spectral, vertically real
      fields.  SINFT is its own inverse.  To find inverse SINFT 
      from 2D-spectra to horizontally spectra, vertically real fields

   forward transform from fld->ffld
   backward transform from ffld->fld

   both fld,ffld have size 0..nhf, 0..mm
   note: sinft(y[],n) expects y to range from 1..n
   =====================================================================
*/
void vsfft2D(fld,nhf,mm, ffld,isgn)
   int nhf,mm,isgn;
   dcomplex **fld,**ffld;
{
   int n,m;
   double norm;
   double *rfld;

   rfld=dvector(1,mm);
   norm = 2.0/((double) mm); /* see comments in header of hfft2D */

   if (isgn==1) {
      for (n=0;n<=nhf;n++) {
         for (m=1;m<=mm;m++) rfld[m] = fld[n][m].r;
         mysinft(rfld,mm);
         for (m=1;m<=mm;m++) ffld[n][m].r = norm*rfld[m];
         ffld[n][0].r = ffld[n][mm].r;

         for (m=1;m<=mm;m++) rfld[m] = fld[n][m].i;
         mysinft(rfld,mm);
         for (m=1;m<=mm;m++) ffld[n][m].i = norm*rfld[m];
         ffld[n][0].r = ffld[n][mm].r;
      }
   } else if (isgn==-1) {
      for (n=0;n<=nhf;n++) {
         for (m=1;m<=mm;m++) rfld[m] = ffld[n][m].r;
         mysinft(rfld,mm);
         for (m=1;m<=mm;m++) fld[n][m].r = rfld[m];
         fld[n][0].r = fld[n][mm].r;

         for (m=1;m<=mm;m++) rfld[m] = ffld[n][m].i;
         mysinft(rfld,mm);
         for (m=1;m<=mm;m++) fld[n][m].i = rfld[m];
         fld[n][0].r = fld[n][mm].r;
      }
   }

   free_dvector(rfld,1,mm);
}

/* =====================================================================

   ****  WARNING WARNING WARNING  - THIS HAS NOT BEEN TESTED !!!!!!!

   SUBROUTINE: VCFFT2D
   Calculate fourier COSINE transform of horizontally spectral, vertically real
      fields.  If sgn=-1 find inverse COSFT from 2D-spectra to horizontally
      spectra, vertically real fields

   forward transform from fld->ffld
   backward transform from ffld->fld

   both fld,ffld have size 0..nhf, 0..mm
   note: cosft(y[],n,isgn) expects y to range from 1..n+1
   =====================================================================
*/
void vcfft2D(fld,nhf,mm, ffld,isgn)
   int nhf,mm,isgn;
   dcomplex **fld,**ffld;
{
   int n,m;
   double norm;
   double *rfld;

   rfld=dvector(0,mm);
   norm = 2.0/((double) mm); /* see comments in header of hfft2D */

   if (isgn==1) {
      for (n=0;n<=nhf;n++) {
         for (m=0;m<=mm;m++) rfld[m] = fld[n][m].r;
         mycosft(rfld-1,mm);
         for (m=0;m<=mm;m++) ffld[n][m].r = norm*rfld[m];

         for (m=0;m<=mm;m++) rfld[m] = fld[n][m].i;
         mycosft(rfld-1,mm);
         for (m=0;m<=mm;m++) ffld[n][m].i = norm*rfld[m];
      }
   } else if (isgn==-1) {
      for (n=0;n<=nhf;n++) {
         for (m=0;m<=mm;m++) rfld[m] = ffld[n][m].r;
         mycosft(rfld-1,mm);
         for (m=0;m<=mm;m++) fld[n][m].r = rfld[m];

         for (m=0;m<=mm;m++) rfld[m] = ffld[n][m].i;
         mycosft(rfld-1,mm);
         for (m=0;m<=mm;m++) fld[n][m].i = rfld[m];
      }
   }

   free_dvector(rfld,0,mm);
}

/* =====================================================================

   SUBROUTINE: SFFT2D    (tested and approved  Aug 7, 93)

   Calculate fourier SINE transform of real 2D field.  sinft is its own
      inverse though a normalization factor 2/nx = 2 (dx/L) should be 
      multiplied to the forward transform so that F represents the
      _average_ and not the sum of components of f.
      
   note: sinft(y[],n) expects y to range from 1..n where y[1]=0.0;
   =====================================================================
*/
void sfft2D(ar,nx,ny,index, far,isgn)
   int nx,ny,index,isgn;
   double **ar,**far;
{
   int i,j;
   double norm;
   double *v;

   if (index==1) {   /* sin transform x component */
     v=dvector(0,nx);
     norm = 2.0/((double) nx); 

     if (isgn==1) { 
        if (nx<=0) myerror("SFFT2D: nx too small for sine transform");
        for (j=0;j<=ny;j++) {
           for (i=0;i<=nx;i++) v[i] = ar[i][j];
           mysinft(v-1,nx);
           v[nx]=0.0;
           for (i=0;i<=nx;i++) far[i][j] = norm*v[i];
        }
     } else if (isgn==-1) {
        if (nx<=0) myerror("SFFT2D: nx too small for sine transform");
        for (j=0;j<=ny;j++) {
           for (i=0;i<=nx;i++) v[i] = far[i][j];
           mysinft(v-1,nx);
           v[nx]=0.0;
           for (i=0;i<=nx;i++) ar[i][j] = v[i];
        }
     }
     free_dvector(v,0,nx);

   } else if (index==2) {  /* sin transform y component */

     v=dvector(0,ny);
     norm = 2.0/((double) ny); 

     if (isgn==1) {
        if (ny<=0) myerror("SFFT2D: ny too small for sine transform");
        for (i=0;i<=nx;i++) {
           for (j=0;j<=ny;j++) v[j] = ar[i][j];
           mysinft(v-1,ny);
           v[ny]=0.0;
           for (j=0;j<=ny;j++) far[i][j] = norm*v[j];
        }
     } else if (isgn==-1) {
        if (ny<=0) myerror("SFFT2D: ny too small for sine transform");
        for (i=0;i<=nx;i++) {
           for (j=0;j<=ny;j++) v[j] = far[i][j];
           mysinft(v-1,ny);
           v[ny]=0.0;
           for (j=0;j<=ny;j++) ar[i][j] = v[j];
        }
     }
     free_dvector(v,0,ny);
   }
}

/* =====================================================================

   SUBROUTINE: CFFT2D    (tested and approved  Aug 8, 93)

   Calculate fourier COSINE transform of real 2D field.  cosft is its own
      inverse though a normalization factor 2/nx = 2 (dx/L) should be 
      multiplied to the forward transform so that F represents the
      _average_ and not the sum of components of f.
      
   note: cosft(y[],n,isgn) expects y to range from 1..n+1, where y[n+1]=y[1]
    
   FEB 13, 2004: Modified so starting indices also provided
   =====================================================================
*/
void cfft2D(ar,ox,nx,oy,ny,index, far,isgn)
   int ox,nx,oy,ny,index,isgn;
   double **ar,**far;
{
   int i,j;
   int n;
   double norm;
   double *v;


   if (index==1) {   /* cos transform x component */
     if (nx<=ox) myerror("CFFT2D: nx too small for cosine transform");
     n=nx-ox;
     v=dvector(0,n);
     norm = 2.0/((double) n); 

     if (isgn==1) { 
       for (j=oy;j<=ny;j++) {
           /* only nx points used in transform */
           for (i=0;i<=n;i++) v[i] = ar[i+ox][j];  
           mycosft(v-1,n);
           v[0] *= 0.5;
           for (i=ox;i<=nx;i++) far[i][j] = norm*v[i-ox];
       }
     } else if (isgn==-1) {
       for (j=oy;j<=ny;j++) {
           /* only nx points is used in transform */
           for (i=0;i<=n;i++) v[i] = far[i+ox][j];
           v[0] *= 2.0;
           mycosft(v-1,n);
           for (i=ox;i<=nx;i++) ar[i][j] = v[i-ox];
       }
     }
     free_dvector(v,0,n);

   } else if (index==2) {  /* cos transform y component */

     if (ny<=oy) myerror("CFFT2D: ny too small for cosine transform");
     n=ny-oy;
     v=dvector(0,n);
     norm = 2.0/((double) n); 

     if (isgn==1) {
       for (i=ox;i<=nx;i++) {
         /* for (j=0;j<=n;j++) v[j] = ar[i][j+oy]; */
         mycosft(ar[i]+oy-1,n);
         /* v[0] *= 0.5; */
         ar[i][oy] *= 0.5;
         for (j=oy;j<=ny;j++) far[i][j] = norm*ar[i][j];
       }
     } else if (isgn==-1) {
       for (i=ox;i<=nx;i++) {
         /* for (j=0;j<=n;j++) v[j] = far[i][j+oy]; */
         far[i][oy] *= 2.0;
         mycosft(far[i]+oy-1,n);
         for (j=oy;j<=ny;j++) ar[i][j] = far[i][j];
       }
     }
     free_dvector(v,0,n);
   }
}

/* ======================================================================
   SUBROUTINE: znfft2D

   This routine is virtually identical to hfft2D except that some smoothing
   near boundaries is omitted.  The routine is called by the 3D simulation
   program.
   ======================================================================
*/
void znfft2D(mat,nn,oy,mm, cfld,isgn)
   int nn,oy,mm,isgn;
   double **mat;
   dcomplex **cfld;
{
   int n,m,n2,n2p1,nhf;
   double norm;
   double *ft;

   nhf = nn>>1;
   norm = 1.0/((double) nn);      /* see comments above */

   /* allocate space */
   ft = dvector(0,nn);

   /* Complex FT values are stored in matrices with components 2n-1 (real)
      and 2n (imaginary).  In frequency space, a[1] is zero frequency
      component (the sum of the real space values) and a[2] is the
      1/2d frequency component (alternating sum of real space values).
   */
   if (isgn==1) {
      for (m=oy;m<=mm;m++) {
         for (n=1;n<=nn;n++) ft[nn-n]=mat[n][m];

         myrealft(ft-1,nn,1);

         cfld[0][m].r = norm*ft[0];
         cfld[nhf][m].r = norm*ft[1];
         cfld[0][m].i = 0.0;
         cfld[nhf][m].i = 0.0;
         for (n=1;n<nhf;n++) {
            n2 = n<<1;
            n2p1 = n2+1;
            cfld[n][m].r = norm*ft[n2];
            cfld[n][m].i = norm*ft[n2p1];
         }
      }
   } else {
      for (m=oy;m<=mm;m++) {  
         ft[0] = cfld[0][m].r;
         ft[1] = cfld[nhf][m].r;
         for (n=1;n<nhf;n++) {
            n2 = n<<1;
            n2p1 = n2+1;
            ft[n2] = cfld[n][m].r;
            ft[n2p1] = cfld[n][m].i;
         }

         /* fft ftvec to physical space */
         myrealft(ft-1,nn,-1);
         for (n=1;n<=nn;n++) mat[n][m] = 2*ft[nn-n];

         /* reproduce first column in last column to show periodicity */
         mat[0][m]=mat[nn][m];
      }
   }

   free_dvector(ft,0,nn);
}

/* ======================================================================
   SUBROUTINE: azfft2D  ("azimuthal fft")

   This routine is similar to znfft2D except that it acts on the azimuthal
   part of a field expressed in polar co-ordinates f(r,theta) returning
   Ff(r,l).
   ======================================================================
*/
void azfft2D(mat,or,nr,nt, cfld,isgn)
   int or,nr,nt,isgn;
   double **mat;
   dcomplex **cfld;
{
   int n,m,m2,m2p1,nthf;
   double norm;
   double *ft;

   nthf = nt>>1;
   norm = 1.0/((double) nt);      /* see comments to hfft2D */

   /* allocate space */
   ft = dvector(0,nt);

   /* Complex FT values are stored in matrices with components 2n-1 (real)
      and 2n (imaginary).  In frequency space, a[1] is zero frequency
      component (the sum of the real space values) and a[2] is the
      1/2d frequency component (alternating sum of real space values).
   */
   if (isgn==1) {
      for (n=or;n<=nr;n++) {
         for (m=1;m<=nt;m++) ft[nt-m]=mat[n][m];

         myrealft(ft-1,nt,1);

         cfld[n][0].r = norm*ft[0];     cfld[n][0].i = 0.0;
         cfld[n][nthf].r = norm*ft[1];  cfld[n][nthf].i = 0.0;
         for (m=1;m<nthf;m++) {
            m2 = m<<1;
            m2p1 = m2+1;
            cfld[n][m].r = norm*ft[m2];
            cfld[n][m].i = norm*ft[m2p1];
         }
      }
   } else {
      for (n=or;n<=nr;n++) {
         ft[0] = cfld[n][0].r;
         ft[1] = cfld[n][nthf].r;
         for (m=1;m<nthf;m++) {
            m2 = m<<1;
            m2p1 = m2+1;
            ft[m2] = cfld[n][m].r;
            ft[m2p1] = cfld[n][m].i;
         }

         /* fft ftvec to physical space */
         myrealft(ft-1,nt,-1);
         for (m=1;m<=nt;m++) mat[n][m] = 2*ft[nt-m];

         /* reproduce first column in last column to show periodicity */
         mat[n][0]=mat[n][nt];
      }
   }

   free_dvector(ft,0,nt);
}


/* ======================================================================
   SUBROUTINE: zhfft2D
   zomplex version of hfft2D
   ======================================================================
*/
void zhfft2D(mat,nn,mm, cfld,isgn)
   int nn,mm,isgn;
   double **mat;
   zomplex **cfld;
{
   int n,m,n2,n2p1,nhf;
   double norm;
   double *ft;

   nhf = nn>>1;
   norm = 1.0/((double) nn);      /* see comments above */

   /* allocate space */
   ft = dvector(0,nn);

   /* Complex FT values are stored in matrices with components 2n-1 (real)
      and 2n (imaginary).  In frequency space, a[1] is zero frequency
      component (the sum of the real space values) and a[2] is the
      1/2d frequency component (alternating sum of real space values).
   */
   if (isgn==1) {
      for (m=0;m<=mm;m++) {
         for (n=1;n<=nn;n++) ft[nn-n]=mat[n][m];

         myrealft(ft-1,nn,1);

         cfld[0][m].re = norm*ft[0];
         cfld[nhf][m].re = norm*ft[1];
         cfld[0][m].im = 0.0;
         cfld[nhf][m].im = 0.0;
         for (n=1;n<nhf;n++) {
            n2 = n<<1;
            n2p1 = n2+1;
            cfld[n][m].re = norm*ft[n2];
            cfld[n][m].im = norm*ft[n2p1];
         }
      }
   } else {
      for (m=0;m<=mm;m++) {  
         ft[0] = cfld[0][m].re;
         ft[1] = cfld[nhf][m].re;
         for (n=1;n<nhf;n++) {
            n2 = n<<1;
            n2p1 = n2+1;
            ft[n2] = cfld[n][m].re;
            ft[n2p1] = cfld[n][m].im;
         }

         /* fft ftvec to physical space */
         myrealft(ft-1,nn,-1);
         for (n=1;n<=nn;n++) mat[n][m] = 2*ft[nn-n];

         /* reproduce first column in last column to show periodicity */
         mat[0][m]=mat[nn][m];
      }
   }

   free_dvector(ft,0,nn);
}

/* ======================================================================
   SUBROUTINE: zchfft2D
   Like zhfft2D but transforming a doubly spectral field
   to produce a vertically spectral, horizontally spatial field
   ======================================================================
*/
void zchfft2D(mat,nn,mm, cfld,isgn)
   int nn,mm,isgn;
   double **mat;
   zomplex **cfld;
{
   int n,m,n2,n2p1,nhf;
   double norm;
   double *ft;

   nhf = nn>>1;
   norm = 1.0/((double) nn);      /* see comments above */

   /* allocate space */
   ft = dvector(0,nn);

   /* Complex FT values are stored in matrices with components 2n-1 (real)
      and 2n (imaginary).  In frequency space, a[1] is zero frequency
      component (the sum of the real space values) and a[2] is the
      1/2d frequency component (alternating sum of real space values).
   */
   if (isgn==1) {
      for (m=0;m<=mm;m++) {
         for (n=1;n<=nn;n++) ft[nn-n]=mat[n][m];

         myrealft(ft-1,nn,1);

         cfld[0][m].re = norm*ft[0];
         cfld[nhf][m].re = norm*ft[1];
         cfld[0][m].im = 0.0;
         cfld[nhf][m].im = 0.0;
         for (n=1;n<nhf;n++) {
            n2 = n<<1;
            n2p1 = n2+1;
            cfld[n][m].re = norm*ft[n2];
            cfld[n][m].im = norm*ft[n2p1];
         }
      }
   } else {
      for (m=0;m<=mm;m++) {  
         ft[0] = cfld[0][m].re;
         ft[1] = cfld[nhf][m].re;
         for (n=1;n<nhf;n++) {
            n2 = n<<1;
            n2p1 = n2+1;
            ft[n2] = cfld[n][m].re;
            ft[n2p1] = cfld[n][m].im;
         }

         /* fft ftvec to physical space */
         myrealft(ft-1,nn,-1);
         for (n=1;n<=nn;n++) mat[n][m] = 2*ft[nn-n];

         /* reproduce first column in last column to show periodicity */
         mat[0][m]=mat[nn][m];
      }
   }

   free_dvector(ft,0,nn);
}

/* =====================================================================
   SUBROUTINE: zVRFFT2D
   Calculate vertical fourier transform of real field.  
      If sgn=-1 find inverse FT from vertically spectral to
      vertically real field

   mat is a nn by mm real matrix ( mat[n][0]==mat[n][mm] )
   cfld is an nn by (mm/2) complex matrix  (cfld[nn][0].im=cfld[nn][mm/2].im=0)
  
   The point mat[0][0] corresponds to x=z=0.
   The element cfld[n][0] corresponds to the zero frequency mode at level n.
   =====================================================================
*/
void zvrfft2D(double **mat,int nn,int mm,zomplex **cfld, int isgn)
{
   int n,m,m2,m2p1,mhf;
   double norm;
   double *ft;

   mhf = mm>>1;
   norm = 1.0/((double) mm);      /* see to HFFT2D comments */

   /* allocate space */
   ft = dvector(0,mm);

   /* Complex FT values are stored in matrices with components 2n-1 (real)
      and 2n (imaginary).  In frequency space, a[1] is zero frequency
      component (the sum of the real space values) and a[2] is the
      1/2d frequency component (alternating sum of real space values).
   */
   if (isgn==1) {
      for (n=0;n<nn;n++) {  /* cfld at nn equals cfld at zero */
         for (m=1;m<=mm;m++) ft[mm-m]=mat[n][m];

         myrealft(ft-1,mm,1);

         cfld[n][0].re = norm*ft[0];
         cfld[n][mhf].re = norm*ft[1];
         cfld[n][0].im = 0.0;
         cfld[n][mhf].im = 0.0;
         for (m=1;m<mhf;m++) {
            m2 = m<<1;
            m2p1 = m2+1;
            cfld[n][m].re = norm*ft[m2];
            cfld[n][m].im = norm*ft[m2p1];
         }
      }
      for (m=0;m<=mhf;m++) cfld[nn][m]=cfld[0][m];
   } else {
      for (n=0;n<=nn;n++) { 
         ft[0] = cfld[n][0].re;
         ft[1] = cfld[n][mhf].re;
         for (m=1;m<mhf;m++) {
            m2 = m<<1;
            m2p1 = m2+1;
            ft[m2] = cfld[n][m].re;
            ft[m2p1] = cfld[n][m].im;
         }

         /* fft ftvec to physical space */
         myrealft(ft-1,mm,-1);
         for (m=1;m<=mm;m++) mat[n][m] = 2*ft[mm-m];

         /* reproduce first column in last column to show periodicity */
         mat[n][0]=mat[n][mm];
      }
   }
   free_dvector(ft,0,mm);
}

/* =====================================================================
   SUBROUTINE: zVFFT2D
   Calculate fourier transform of horizontally spectral, vertically real
      fields.  If sgn=-1 find inverse FT from 2D-spectra to horizontally
      spectral, vertically real fields

   forward transform from fld->ffld
   backward transform from ffld->fld

   The element fld[0][m] corresponds to the zero horiz. freq at level m.
   The element ffld[0][mhf] corresponds to the zero frequency mode

   both fld,ffld have size 0..nhf, 0..mm
   =====================================================================
*/
void zvfft2D(fld,nhf,mm, ffld,isgn)
   int nhf,mm,isgn;
   zomplex **fld,**ffld;
{
   int mhf,n,m,m2,mpmhf,mmmhf,mmm;
   double norm;
   double *rfld;

   rfld=dvector(0,2*mm);

   mhf = mm>>1;
   norm = 1.0/((double) mm);     /* see comments in header of hfft2D */

   if (isgn==1) {
      for (n=0;n<=nhf;n++) {
         for (m=0;m<mm;m++) {
            m2 = m<<1;
            mmm=mm-m;
            rfld[m2] = fld[n][mmm].re;
            rfld[m2+1] = fld[n][mmm].im;
         }

         myfour1(rfld-1,mm,1);

         for (m=0;m<mhf;m++) {
            m2=m<<1;
            mpmhf=mhf+m;
            ffld[n][mpmhf].re = norm*rfld[m2];
            ffld[n][mpmhf].im = norm*rfld[m2+1];
         }
         for (m=mhf;m<mm;m++) {
            m2=m<<1;
            mmmhf= m-mhf;
            ffld[n][mmmhf].re = norm*rfld[m2];
            ffld[n][mmmhf].im = norm*rfld[m2+1];
         }
         /* enforce zero values of im parts for real vectors */
         if (n==0 || n==nhf) ffld[n][mhf].im = ffld[n][0].im = 0.0;

         ffld[n][mm]=ffld[n][0];
      }
   } else if (isgn==-1) {
      for (n=0;n<=nhf;n++) {
         for (m=0;m<mhf;m++) {
            m2=m<<1;
            mpmhf=m+mhf;
            rfld[m2] = ffld[n][mpmhf].re;
            rfld[m2+1] = ffld[n][mpmhf].im;
         }
         for (m=mhf;m<mm;m++) {
            m2=m<<1;
            mmmhf=m-mhf;
            rfld[m2] = ffld[n][mmmhf].re;
            rfld[m2+1] = ffld[n][mmmhf].im;
         }

         myfour1(rfld-1,mm,-1);

         for (m=0;m<mm;m++) {
            m2=m<<1;
            mmm=mm-m;
            fld[n][mmm].re = rfld[m2];
            fld[n][mmm].im = rfld[m2+1];
         }
         /* enforce zero values of im parts for real vectors */
         if (n==0 || n==nhf) for(m=0;m<=mm;m++) fld[n][m].im = 0.0;

         fld[n][0] = fld[n][mm];
      }
   }

   free_dvector(rfld,0,2*mm);
}

/* =====================================================================

   ****  WARNING WARNING WARNING  - THIS HAS NOT BEEN TESTED !!!!!!!

   SUBROUTINE: zVSFFT2D
   Calculate fourier SINE transform of horizontally spectral, vertically real
      fields.  SINFT is its own inverse.  To find inverse SINFT 
      from 2D-spectra to horizontally spectra, vertically real fields

   forward transform from fld->ffld
   backward transform from ffld->fld

   both fld,ffld have size 0..nhf, 0..mm
   note: sinft(y[],n) expects y to range from 1..n
   =====================================================================
*/
void zvsfft2D(fld,nhf,mm, ffld,isgn)
   int nhf,mm,isgn;
   zomplex **fld,**ffld;
{
   int n,m;
   double norm;
   double *rfld;

   rfld=dvector(1,mm);
   norm = 2.0/((double) mm); /* see comments in header of hfft2D */

   if (isgn==1) {
      for (n=0;n<=nhf;n++) {
         for (m=1;m<=mm;m++) rfld[m] = fld[n][m].re;
         mysinft(rfld,mm);
         for (m=1;m<=mm;m++) ffld[n][m].re = norm*rfld[m];
         ffld[n][0].re = ffld[n][mm].re;

         for (m=1;m<=mm;m++) rfld[m] = fld[n][m].im;
         mysinft(rfld,mm);
         for (m=1;m<=mm;m++) ffld[n][m].im = norm*rfld[m];
         ffld[n][0].re = ffld[n][mm].re;
      }
   } else if (isgn==-1) {
      for (n=0;n<=nhf;n++) {
         for (m=1;m<=mm;m++) rfld[m] = ffld[n][m].re;
         mysinft(rfld,mm);
         for (m=1;m<=mm;m++) fld[n][m].re = rfld[m];
         fld[n][0].re = fld[n][mm].re;

         for (m=1;m<=mm;m++) rfld[m] = ffld[n][m].im;
         mysinft(rfld,mm);
         for (m=1;m<=mm;m++) fld[n][m].im = rfld[m];
         fld[n][0].re = fld[n][mm].re;
      }
   }

   free_dvector(rfld,1,mm);
}

/* =====================================================================

   ****  WARNING WARNING WARNING  - THIS HAS NOT BEEN TESTED !!!!!!!

   SUBROUTINE: zVCFFT2D
   Calculate fourier COSINE transform of horizontally spectral, vertically real
      fields.  If sgn=-1 find inverse COSFT from 2D-spectra to horizontally
      spectra, vertically real fields

   forward transform from fld->ffld
   backward transform from ffld->fld

   both fld,ffld have size 0..nhf, 0..mm
   note: cosft(y[],n,isgn) expects y to range from 1..n+1
   =====================================================================
*/
void zvcfft2D(fld,nhf,mm, ffld,isgn)
   int nhf,mm,isgn;
   zomplex **fld,**ffld;
{
   int n,m;
   double norm;
   double *rfld;

   rfld=dvector(0,mm);
   norm = 2.0/((double) mm); /* see comments in header of hfft2D */

   if (isgn==1) {
      for (n=0;n<=nhf;n++) {
         for (m=0;m<=mm;m++) rfld[m] = fld[n][m].re;
         mycosft(rfld-1,mm);
         for (m=0;m<=mm;m++) ffld[n][m].re = norm*rfld[m];

         for (m=0;m<=mm;m++) rfld[m] = fld[n][m].im;
         mycosft(rfld-1,mm);
         for (m=0;m<=mm;m++) ffld[n][m].im = norm*rfld[m];
      }
   } else if (isgn==-1) {
      for (n=0;n<=nhf;n++) {
         for (m=0;m<=mm;m++) rfld[m] = ffld[n][m].re;
         mycosft(rfld-1,mm);
         for (m=0;m<=mm;m++) fld[n][m].re = rfld[m];

         for (m=0;m<=mm;m++) rfld[m] = ffld[n][m].im;
         mycosft(rfld-1,mm);
         for (m=0;m<=mm;m++) fld[n][m].im = rfld[m];
      }
   }

   free_dvector(rfld,0,mm);
}

/* ======================================================================
   SUBROUTINE: zznfft2D

   This routine is virtually identical to hfft2D except that some smoothing
   near boundaries is omitted.  The routine is called by the 3D simulation
   program.
   ======================================================================
*/
void zznfft2D(mat,nn,oy,mm, cfld,isgn)
   int nn,oy,mm,isgn;
   double **mat;
   zomplex **cfld;
{
   int n,m,n2,n2p1,nhf;
   double norm;
   double *ft;

   nhf = nn>>1;
   norm = 1.0/((double) nn);      /* see comments above */

   /* allocate space */
   ft = dvector(0,nn);

   /* Complex FT values are stored in matrices with components 2n-1 (real)
      and 2n (imaginary).  In frequency space, a[1] is zero frequency
      component (the sum of the real space values) and a[2] is the
      1/2d frequency component (alternating sum of real space values).
   */
   if (isgn==1) {
      for (m=oy;m<=mm;m++) {
         for (n=1;n<=nn;n++) ft[nn-n]=mat[n][m];

         myrealft(ft-1,nn,1);

         cfld[0][m].re = norm*ft[0];
         cfld[nhf][m].re = norm*ft[1];
         cfld[0][m].im = 0.0;
         cfld[nhf][m].im = 0.0;
         for (n=1;n<nhf;n++) {
            n2 = n<<1;
            n2p1 = n2+1;
            cfld[n][m].re = norm*ft[n2];
            cfld[n][m].im = norm*ft[n2p1];
         }
      }
   } else {
      for (m=oy;m<=mm;m++) {  
         ft[0] = cfld[0][m].re;
         ft[1] = cfld[nhf][m].re;
         for (n=1;n<nhf;n++) {
            n2 = n<<1;
            n2p1 = n2+1;
            ft[n2] = cfld[n][m].re;
            ft[n2p1] = cfld[n][m].im;
         }

         /* fft ftvec to physical space */
         myrealft(ft-1,nn,-1);
         for (n=1;n<=nn;n++) mat[n][m] = 2*ft[nn-n];

         /* reproduce first column in last column to show periodicity */
         mat[0][m]=mat[nn][m];
      }
   }

   free_dvector(ft,0,nn);
}

/* ======================================================================
   SUBROUTINE: zazfft2D  ("azimuthal fft")

   This routine is similar to znfft2D except that it acts on the azimuthal
   part of a field expressed in polar co-ordinates f(r,theta) returning
   Ff(r,l).
   ======================================================================
*/
void zazfft2D(mat,or,nr,nt, cfld,isgn)
   int or,nr,nt,isgn;
   double **mat;
   zomplex **cfld;
{
   int n,m,m2,m2p1,nthf;
   double norm;
   double *ft;

   nthf = nt>>1;
   norm = 1.0/((double) nt);      /* see comments to hfft2D */

   /* allocate space */
   ft = dvector(0,nt);

   /* Complex FT values are stored in matrices with components 2n-1 (real)
      and 2n (imaginary).  In frequency space, a[1] is zero frequency
      component (the sum of the real space values) and a[2] is the
      1/2d frequency component (alternating sum of real space values).
   */
   if (isgn==1) {
      for (n=or;n<=nr;n++) {
         for (m=1;m<=nt;m++) ft[nt-m]=mat[n][m];

         myrealft(ft-1,nt,1);

         cfld[n][0].re = norm*ft[0];     cfld[n][0].im = 0.0;
         cfld[n][nthf].re = norm*ft[1];  cfld[n][nthf].im = 0.0;
         for (m=1;m<nthf;m++) {
            m2 = m<<1;
            m2p1 = m2+1;
            cfld[n][m].re = norm*ft[m2];
            cfld[n][m].im = norm*ft[m2p1];
         }
      }
   } else {
      for (n=or;n<=nr;n++) {
         ft[0] = cfld[n][0].re;
         ft[1] = cfld[n][nthf].re;
         for (m=1;m<nthf;m++) {
            m2 = m<<1;
            m2p1 = m2+1;
            ft[m2] = cfld[n][m].re;
            ft[m2p1] = cfld[n][m].im;
         }

         /* fft ftvec to physical space */
         myrealft(ft-1,nt,-1);
         for (m=1;m<=nt;m++) mat[n][m] = 2*ft[nt-m];

         /* reproduce first column in last column to show periodicity */
         mat[n][0]=mat[n][nt];
      }
   }

   free_dvector(ft,0,nt);
}
