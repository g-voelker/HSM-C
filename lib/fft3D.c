/* =========================================================================
   Routines to take fourier transforms of 3-D fields

   znfft3D : take fourier/inverse fourier transform in zonal direction
             assume periodic and uses exp(I n alpha x) as basis
   mdfft3D : take fourier transform in meridional direction
   vfft3D  : take fourier transform in vertical direction
   =========================================================================
*/

#include <math.h>
#include "complex.h"
#include "alloc_space.h"
#include "myerror.h"
#include "myfft.h"

/* ======================================================================
   SUBROUTINE: znfft3D
   Calculate fourier transform of real field to zonal spectral, meridionally 
      and vertically real field.  If sgn=-1 find inverse FT from horizontally
      spectral, meridionally and vertically real field to real field.  Note - 
      for a real field the fourier transform has negative frequency components 
      which are complex conjugates of the positive frequency parts.  Only
      the positive frequency components are returned.

   f is an nx,ny,nz real matrix (mat[0][j][k]==mat[nx][j][k])
   F is an (nnx=nx/2),ny,nz complex matrix (F[0][j][k].i=F[nnx][j][k].i=0)
  
   The point f[0][0][0] corresponds to x=y=z=0.
   The element F[0][j][k] corresponds to the zero frequency zonal mode.

   According to Numerical Recipes, fourier transforms of nx data to
      nnx complex frequency components are not normalized but inverse
      transform is normalized by  nnx.  In this case the zero frequency
      term will correspond to the sum of all real data (NOT THE AVERAGE!).
      Well, this sucks.  It makes much more sense to take the average by
      dividing the forward transform by nx.  So that's what I do.  Nyah. 
   ======================================================================
*/
void znfft3D(f,nx,oy,ny,oz,nz, F,isgn)
   int nx,oy,ny,oz,nz,isgn;
   double ***f;
   dcomplex ***F;
{
   int n,j,k,n2,nnx;
   double norm;
   double *ft;

   nnx = nx>>1;
   norm = 1.0/((double) nx);      /* see comments above */

   /* allocate space */
   ft = dvector(0,nx);

   /* Complex FT values are stored in matrices with components 2n-1 (real)
      and 2n (imaginary).  In frequency space, a[1] is zero frequency
      component (the sum of the real space values) and a[2] is the
      1/2d frequency component (alternating sum of real space values).
   */
   if (isgn==1) {
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) {
         for (n=1;n<=nx;n++) ft[nx-n]=f[n][j][k];

         myrealft(ft-1,nx,1);
   
         F[0][j][k].r = norm*ft[0];
         F[nnx][j][k].r = norm*ft[1];
         F[0][j][k].i = 0.0;
         F[nnx][j][k].i = 0.0;
         for (n=1;n<nnx;n++) {
            n2 = n<<1;
            F[n][j][k].r = norm*ft[n2];
            F[n][j][k].i = norm*ft[n2+1];
         }
      }
   } else {
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) {
         ft[0] = F[0][j][k].r;
         ft[1] = F[nnx][j][k].r;
         for (n=1;n<nnx;n++) {
            n2 = n<<1;
            ft[n2] = F[n][j][k].r;
            ft[n2+1] = F[n][j][k].i;
         }

         /* fft ftvec to physical space */
         myrealft(ft-1,nx,-1);
         for (n=1;n<=nx;n++) f[n][j][k] = 2*ft[nx-n];

         /* reproduce first column in last column to show periodicity */
         f[0][j][k]=f[nx][j][k];
      }
   }

   free_dvector(ft,0,nx);
}

/* =====================================================================
   SUBROUTINE: MDFFT3D
   Calculate fourier transform of zonally spectral, meridionally real
      fields.  If sgn=-1 find inverse FT from 3D-spectra to zonally
      spectral, meridionally real fields

      isgn=1 :  F -> cF        isgn=-1 :  cF -> F

   The element F[0][j][k] corresponds to the zero horiz. freq at latitude j.
   The element cF[0][nny][k] corresponds to the zero frequency mode, nny=ny/2

   both F and cF have size 0..nnx, 0..ny, 0..nz
   =====================================================================
*/
void mdfft3D(F,nnx,ny,oz,nz, cF,isgn)
   int nnx,ny,oz,nz,isgn;
   dcomplex ***F,***cF;
{
   int n,j,k,j2,tj,nny,ny2;
   double norm;
   double *rfld;

   if (ny<=0) myerror("MDFFT3D: ny too small for fourier transform");

   ny2 = ny<<1;
   rfld=dvector(0,ny2-1);

   nny = ny>>1;
   norm = 1.0/((double) ny);     /* see comments in header of hfft2D */

   if (isgn==1) {
      for (n=0;n<=nnx;n++) 
      for (k=oz;k<=nz;k++) {
         for (j=0;j<ny;j++) {
            j2 = j<<1;
            tj=ny-j;
            rfld[j2] = F[n][tj][k].r;
            rfld[j2+1] = F[n][tj][k].i;
         }

         myfour1(rfld-1,ny,1);

         for (j=0;j<nny;j++) {
            j2=j<<1;
            tj=nny+j;
            cF[n][tj][k].r = norm*rfld[j2];
            cF[n][tj][k].i = norm*rfld[j2+1];
         }
         for (j=nny;j<ny;j++) {
            j2=j<<1;
            tj= j-nny;
            cF[n][tj][k].r = norm*rfld[j2];
            cF[n][tj][k].i = norm*rfld[j2+1];
         }
         /* enforce zero values of im parts for real vectors */
         if (n==0 || n==nnx) cF[n][nny][k].i = cF[n][0][k].i = 0.0;

         cF[n][ny][k]=cF[n][0][k];
      }
   } else if (isgn==-1) {
      for (n=0;n<=nnx;n++) 
      for (k=oz;k<=nz;k++) {
         for (j=0;j<nny;j++) {
            j2=j<<1;
            tj=j+nny;
            rfld[j2] = cF[n][tj][k].r;
            rfld[j2+1] = cF[n][tj][k].i;
         }
         for (j=nny;j<ny;j++) {
            j2=j<<1;
            tj=j-nny;
            rfld[j2] = cF[n][tj][k].r;
            rfld[j2+1] = cF[n][tj][k].i;
         }

         myfour1(rfld-1,ny,-1);

         for (j=0;j<ny;j++) {
            j2=j<<1;
            tj=ny-j;
            F[n][tj][k].r = rfld[j2];
            F[n][tj][k].i = rfld[j2+1];
         }
         /* enforce zero values of im parts for real vectors */
         if (n==0 || n==nnx) for(j=0;j<=ny;j++) F[n][j][k].i = 0.0;

         F[n][0][k] = F[n][ny][k];
      }
   }

   free_dvector(rfld,0,ny2);
}

/* =====================================================================
   SUBROUTINE: VFFT3D
   Calculate fourier cos transform of horizontally spectral, vertically real
      fields.  If sgn=-1 find inverse FT from 3D-spectra to horizontally
      spectra, vertically real fields

      isgn=1 :  F -> cF        isgn=-1 :  cF -> F

   The element F[0][j][k] corresponds to the zero horiz. freq at level k.
   The element cF[0][j][nnz] corresponds to the zero frequency mode, nnz=nz/2

   both F and cF have size 0..nnx, 0..ny, 0..nz
   =====================================================================
*/
void vfft3D(F,nnx,ny,nz, cF,isgn)
   int nnx,ny,nz,isgn;
   dcomplex ***F,***cF;
{
   int n,j,k,k2,tk,nnz;
   double norm;
   double *rfld;

   if (nz<=0) myerror("VFFT3D: nz too small for cosine transform");

   rfld=dvector(0,nz<<1);

   nnz = nz>>1;
   norm = 1.0/((double) nz);     /* see comments in header of hfft2D */

   if (isgn==1) {
      for (n=0;n<=nnx;n++) 
      for (j=0;j<=ny;j++) {
         for (k=0;k<nz;k++) {
            k2 = k<<1;
			tk=nz-k;
            rfld[k2] = F[n][j][tk].r;
            rfld[k2+1] = F[n][j][tk].i;
         }

         myfour1(rfld-1,nz,1);

         for (k=0;k<nnz;k++) {
            k2=k<<1;
            tk=nnz+k;
            cF[n][j][tk].r = norm*rfld[k2];
            cF[n][j][tk].i = norm*rfld[k2+1];
         }
         for (k=nnz;k<nz;k++) {
            k2=k<<1;
            tk= k-nnz;
            cF[n][j][tk].r = norm*rfld[k2];
            cF[n][j][tk].i = norm*rfld[k2+1];
         }
         /* enforce zero values of im parts for real vectors */
         if (n==0 || n==nnx) cF[n][j][nnz].i = cF[n][j][0].i = 0.0;

         cF[n][j][nz]=cF[n][j][0];
      }
   } else if (isgn==-1) {
      for (n=0;n<=nnx;n++)
      for (j=0;j<=ny;j++) {
         for (k=0;k<nnz;k++) {
            k2=k<<1;
            tk=k+nnz;
            rfld[k2] = cF[n][j][tk].r;
            rfld[k2+1] = cF[n][j][tk].i;
         }
         for (k=nnz;k<nz;k++) {
            k2=k<<1;
            tk=k-nnz;
            rfld[k2] = cF[n][j][tk].r;
            rfld[k2+1] = cF[n][j][tk].i;
         }

         myfour1(rfld-1,nz,-1);

         for (k=0;k<nz;k++) {
            k2=k<<1;
            tk=nz-k;
            F[n][j][tk].r = rfld[k2];
            F[n][j][tk].i = rfld[k2+1];
         }
         /* enforce zero values of im parts for real vectors */
         if (n==0 || n==nnx) for(k=0;k<=nz;k++) F[n][j][k].i = 0.0;

         F[n][j][0] = F[n][j][nz];
      }
   }

   free_dvector(rfld,0,nz<<1);
}


/* ======================================================================
   SUBROUTINE: zznfft3D
   zomplex version of znfft3D
   ======================================================================
*/
void zznfft3D(f,nx,oy,ny,oz,nz, F,isgn)
   int nx,oy,ny,oz,nz,isgn;
   double ***f;
   zomplex ***F;
{
   int n,j,k,n2,nnx;
   double norm;
   double *ft;

   nnx = nx>>1;
   norm = 1.0/((double) nx);      /* see comments above */

   /* allocate space */
   ft = dvector(0,nx);

   /* Complex FT values are stored in matrices with components 2n-1 (real)
      and 2n (imaginary).  In frequency space, a[1] is zero frequency
      component (the sum of the real space values) and a[2] is the
      1/2d frequency component (alternating sum of real space values).
   */
   if (isgn==1) {
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) {
         for (n=1;n<=nx;n++) ft[nx-n]=f[n][j][k];

         myrealft(ft-1,nx,1);
   
         F[0][j][k].re = norm*ft[0];
         F[nnx][j][k].re = norm*ft[1];
         F[0][j][k].im = 0.0;
         F[nnx][j][k].im = 0.0;
         for (n=1;n<nnx;n++) {
            n2 = n<<1;
            F[n][j][k].re = norm*ft[n2];
            F[n][j][k].im = norm*ft[n2+1];
         }
      }
   } else {
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) {
         ft[0] = F[0][j][k].re;
         ft[1] = F[nnx][j][k].re;
         for (n=1;n<nnx;n++) {
            n2 = n<<1;
            ft[n2] = F[n][j][k].re;
            ft[n2+1] = F[n][j][k].im;
         }

         /* fft ftvec to physical space */
         myrealft(ft-1,nx,-1);
         for (n=1;n<=nx;n++) f[n][j][k] = 2*ft[nx-n];

         /* reproduce first column in last column to show periodicity */
         f[0][j][k]=f[nx][j][k];
      }
   }

   free_dvector(ft,0,nx);
}

/* =====================================================================
   SUBROUTINE: MDFFT3D
   Calculate fourier transform of zonally spectral, meridionally real
      fields.  If sgn=-1 find inverse FT from 3D-spectra to zonally
      spectral, meridionally real fields

      isgn=1 :  F -> cF        isgn=-1 :  cF -> F

   The element F[0][j][k] corresponds to the zero horiz. freq at latitude j.
   The element cF[0][nny][k] corresponds to the zero frequency mode, nny=ny/2

   both F and cF have size 0..nnx, 0..ny, 0..nz
   =====================================================================
*/
void zmdfft3D(F,nnx,ny,oz,nz, cF,isgn)
   int nnx,ny,oz,nz,isgn;
   zomplex ***F,***cF;
{
   int n,j,k,j2,tj,nny,ny2;
   double norm;
   double *rfld;

   if (ny<=0) myerror("MDFFT3D: ny too small for fourier transform");

   ny2 = ny<<1;
   rfld=dvector(0,ny2-1);

   nny = ny>>1;
   norm = 1.0/((double) ny);     /* see comments in header of hfft2D */

   if (isgn==1) {
      for (n=0;n<=nnx;n++) 
      for (k=oz;k<=nz;k++) {
         for (j=0;j<ny;j++) {
            j2 = j<<1;
            tj=ny-j;
            rfld[j2] = F[n][tj][k].re;
            rfld[j2+1] = F[n][tj][k].im;
         }

         myfour1(rfld-1,ny,1);

         for (j=0;j<nny;j++) {
            j2=j<<1;
            tj=nny+j;
            cF[n][tj][k].re = norm*rfld[j2];
            cF[n][tj][k].im = norm*rfld[j2+1];
         }
         for (j=nny;j<ny;j++) {
            j2=j<<1;
            tj= j-nny;
            cF[n][tj][k].re = norm*rfld[j2];
            cF[n][tj][k].im = norm*rfld[j2+1];
         }
         /* enforce zero values of im parts for real vectors */
         if (n==0 || n==nnx) cF[n][nny][k].im = cF[n][0][k].im = 0.0;

         cF[n][ny][k]=cF[n][0][k];
      }
   } else if (isgn==-1) {
      for (n=0;n<=nnx;n++) 
      for (k=oz;k<=nz;k++) {
         for (j=0;j<nny;j++) {
            j2=j<<1;
            tj=j+nny;
            rfld[j2] = cF[n][tj][k].re;
            rfld[j2+1] = cF[n][tj][k].im;
         }
         for (j=nny;j<ny;j++) {
            j2=j<<1;
            tj=j-nny;
            rfld[j2] = cF[n][tj][k].re;
            rfld[j2+1] = cF[n][tj][k].im;
         }

         myfour1(rfld-1,ny,-1);

         for (j=0;j<ny;j++) {
            j2=j<<1;
            tj=ny-j;
            F[n][tj][k].re = rfld[j2];
            F[n][tj][k].im = rfld[j2+1];
         }
         /* enforce zero values of im parts for real vectors */
         if (n==0 || n==nnx) for(j=0;j<=ny;j++) F[n][j][k].im = 0.0;

         F[n][0][k] = F[n][ny][k];
      }
   }

   free_dvector(rfld,0,ny2);
}

/* =====================================================================
   SUBROUTINE: VFFT3D
   Calculate fourier cos transform of horizontally spectral, vertically real
      fields.  If sgn=-1 find inverse FT from 3D-spectra to horizontally
      spectra, vertically real fields

      isgn=1 :  F -> cF        isgn=-1 :  cF -> F

   The element F[0][j][k] corresponds to the zero horiz. freq at level k.
   The element cF[0][j][nnz] corresponds to the zero frequency mode, nnz=nz/2

   both F and cF have size 0..nnx, 0..ny, 0..nz
   =====================================================================
*/
void zvfft3D(F,nnx,ny,nz, cF,isgn)
   int nnx,ny,nz,isgn;
   zomplex ***F,***cF;
{
   int n,j,k,k2,tk,nnz;
   double norm;
   double *rfld;

   if (nz<=0) myerror("VFFT3D: nz too small for cosine transform");

   rfld=dvector(0,nz<<1);

   nnz = nz>>1;
   norm = 1.0/((double) nz);     /* see comments in header of hfft2D */

   if (isgn==1) {
      for (n=0;n<=nnx;n++) 
      for (j=0;j<=ny;j++) {
         for (k=0;k<nz;k++) {
            k2 = k<<1;
			tk=nz-k;
            rfld[k2] = F[n][j][tk].re;
            rfld[k2+1] = F[n][j][tk].im;
         }

         myfour1(rfld-1,nz,1);

         for (k=0;k<nnz;k++) {
            k2=k<<1;
            tk=nnz+k;
            cF[n][j][tk].re = norm*rfld[k2];
            cF[n][j][tk].im = norm*rfld[k2+1];
         }
         for (k=nnz;k<nz;k++) {
            k2=k<<1;
            tk= k-nnz;
            cF[n][j][tk].re = norm*rfld[k2];
            cF[n][j][tk].im = norm*rfld[k2+1];
         }
         /* enforce zero values of im parts for real vectors */
         if (n==0 || n==nnx) cF[n][j][nnz].im = cF[n][j][0].im = 0.0;

         cF[n][j][nz]=cF[n][j][0];
      }
   } else if (isgn==-1) {
      for (n=0;n<=nnx;n++)
      for (j=0;j<=ny;j++) {
         for (k=0;k<nnz;k++) {
            k2=k<<1;
            tk=k+nnz;
            rfld[k2] = cF[n][j][tk].re;
            rfld[k2+1] = cF[n][j][tk].im;
         }
         for (k=nnz;k<nz;k++) {
            k2=k<<1;
            tk=k-nnz;
            rfld[k2] = cF[n][j][tk].re;
            rfld[k2+1] = cF[n][j][tk].im;
         }

         myfour1(rfld-1,nz,-1);

         for (k=0;k<nz;k++) {
            k2=k<<1;
            tk=nz-k;
            F[n][j][tk].re = rfld[k2];
            F[n][j][tk].im = rfld[k2+1];
         }
         /* enforce zero values of im parts for real vectors */
         if (n==0 || n==nnx) for(k=0;k<=nz;k++) F[n][j][k].im = 0.0;

         F[n][j][0] = F[n][j][nz];
      }
   }

   free_dvector(rfld,0,nz<<1);
}
