/* EISPACK routines to find eigenvalues and eigenvectors of a square 
   complex matrix */
#include <math.h>
#include "complex.h"
#include "alloc_space.h"
#include "constants.h"
#include "macros.h"

/* find sqrt(a**2+b**2) without overflow or destructive underflow */
float pythag(float a,float b)
{
   float p,r,s,t,u,tmp;

   p = max(dabs(a),dabs(b));
   if (p!=0.0) {
      tmp = min(dabs(a),dabs(b))/p;
      r = sqr(tmp);
      for (;;) {
         t = 4.0+r;
         if (t==4.0) break;
         s = r/t;
         u = 1.0+2.0*s;
         p *= u;
         r *= sqr(s/u);
         }
      }
   return(p);
}

cdiv( float ar, float ai, float br, float bi, float *cr, float *ci)
{
   float s, ars, ais, brs, bis;

   s = dabs(br) + dabs(bi);
   ars = ar/s;
   ais = ai/s;
   brs = br/s;
   bis = bi/s;

   s = sqr(brs) + sqr(bis);
   *cr = (ars*brs + ais*bis)/s;
   *ci = (ais*brs - ars*bis)/s;
}

csqrt( float xr, float xi, float *yr, float *yi )
{
   float s,tr,ti;

   tr = xr;
   ti = xi;
   s = sqrt( 0.5*(pythag(tr,ti)+dabs(tr)) );
  
   if (tr>=0.0) *yr = s;
   if (ti<0.0) s = -s;
   if (tr<=0.0) *yi=s;
   if (tr<0.0) *yr = 0.5*ti/(*yi);
   if (tr>0.0) *yi = 0.5*ti/(*yr);
}

mswap( float **a, int i1, int j1, int i2, int j2 )
{
   float tmp;

   tmp = a[i1][j1];
   a[i1][j1] = a[i2][j2];
   a[i2][j2] = tmp;
}


comqr2(int n,int low,int high,
       float *ortr,float *orti,float **hr,float **hi,
       float *wr, float *wi,float **zr,float **zi)
/*
      this subroutine finds the eigenvalues and eigenvectors
      of a complex upper hessenberg matrix by the qr
      method.  the eigenvectors of a complex general matrix
      can also be found if  corth  has been used to reduce
      this general matrix to hessenberg form.
 
      on input
         n is the order of the matrix.
         low and high are integers determined by the balancing
           subroutine  cbal.  if  cbal  has not been used,
           set low=1, high=n.
         ortr and orti contain information about the unitary trans-
           formations used in the reduction by  corth, if performed.
           only elements low through high are used.  if the eigenvectors
           of the hessenberg matrix are desired, set ortr[j] and
           orti[j] to 0.0 for these elements.
         hr and hi contain the real and imaginary parts,
           respectively, of the complex upper hessenberg matrix.
           their lower triangles below the subdiagonal contain further
           information about the transformations which were used in the
           reduction by  corth, if performed.  if the eigenvectors of
           the hessenberg matrix are desired, these elements may be
           arbitrary.
 
      on output
         ortr, orti, and the upper hessenberg portions of hr and hi
           have been destroyed.
         wr and wi contain the real and imaginary parts,
           respectively, of the eigenvalues.
         zr and zi contain the real and imaginary parts,
           respectively, of the eigenvectors.  the eigenvectors
           are unnormalized.  if an error exit is made, none of
           the eigenvectors has been found.
*/ 
{
   int i,j,k,l,m,en,ll,itn,its,enm1,iend;
   float si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr;
   float norm,tst1,tst2;
   
   /* initialize eigenvector matrix */
   for (j=1;j<=n;j++) {
      for (i=1;i<=n;i++) {
         zr[i][j] = 0.0;
         zi[i][j] = 0.0;
         }
      zr[j][j] = 1.0;
      }
   
   /* form matrix of accumulated transformations from info left by corth */
   iend = high-low-1;
   if (iend>=0) {
      if (iend!=0) 
         for (i=high-1;i>=low+1;i--) {
            if (ortr[i]!=0.0||orti[i]!=0.0)
               if (hr[i][i-1]!=0.0||hi[i][i-1]!=0.0) {
                  /* norm below is negative of h formed in corth */
                  norm = hr[i][i-1]*ortr[i]+hi[i][i-1]*orti[i];
   
                  for (k=i+1;k<=high;k++) {
                     ortr[k] = hr[k][i-1];
                     orti[k] = hi[k][i-1];
                     }
   
                  for (j=i;j<=high;j++) {
                     sr = 0.0;
                     si = 0.0;
   
                     for (k=i;k<=high;k++) {
                        sr = sr+ortr[k]*zr[k][j]+orti[k]*zi[k][j];
                        si = si+ortr[k]*zi[k][j]-orti[k]*zr[k][j];
                        }
   
                     sr = sr/norm;
                     si = si/norm;
   
                     for (k=i;k<=high;k++) {
                        zr[k][j] = zr[k][j]+sr*ortr[k]-si*orti[k];
                        zi[k][j] = zi[k][j]+sr*orti[k]+si*ortr[k];
                        }
                     }
                  }
            }
      /* create real subdiagonal elements */
      l = low+1;
   
      for (i=l;i<=high;i++) {
         ll = min(i+1,high);
         if (hi[i][i-1]!=0.0) {
            norm = pythag(hr[i][i-1],hi[i][i-1]);
            yr = hr[i][i-1]/norm;
            yi = hi[i][i-1]/norm;
            hr[i][i-1] = norm;
            hi[i][i-1] = 0.0;
   
            for (j=i;j<=n;j++) {
               si = yr*hi[i][j]-yi*hr[i][j];
               hr[i][j] = yr*hr[i][j]+yi*hi[i][j];
               hi[i][j] = si;
               }
   
            for (j=1;j<=ll;j++) {
               si = yr*hi[j][i]+yi*hr[j][i];
               hr[j][i] = yr*hr[j][i]-yi*hi[j][i];
               hi[j][i] = si;
               }
   
            for (j=low;j<=high;j++) {
               si = yr*zi[j][i]+yi*zr[j][i];
               zr[j][i] = yr*zr[j][i]-yi*zi[j][i];
               zi[j][i] = si;
               }
            }
         }
      }
   /* store roots isolated by cbal */
   for (i=1;i<=n;i++) 
      if (i<low||i>high) {
         wr[i] = hr[i][i];
         wi[i] = hi[i][i];
         }
   
   en = high;
   tr = 0.0;
   ti = 0.0;
   itn = 30*n;
   for(;;) {
      /* search for next eigenvalue */
      if (en<low) break;
      its = 0;
      enm1 = en-1;
      for(;;) {
         /* look for a small sub-diagonal element for l=en by -1 until low do */
         for (l=en;l>=low;l--) {
            if (l==low) break;
            tst1 = dabs(hr[l-1][l-1])+dabs(hi[l-1][l-1])+
                   dabs(hr[l][l])+dabs(hi[l][l]);
            tst2 = tst1+dabs(hr[l][l-1]);
            if (tst2==tst1) break;
            }
   
         /* form shift */
         if (l==en || itn==0) break;
         if (its==10||its==20) {
            /* form exceptional shift */
            sr = dabs(hr[en][enm1])+dabs(hr[enm1][en-2]);
            si = 0.0;
         } else {
            sr = hr[en][en];
            si = hi[en][en];
            xr = hr[enm1][en]*hr[en][enm1];
            xi = hi[enm1][en]*hr[en][enm1];
            if (xr!=0.0||xi!=0.0) {
               yr = (hr[enm1][enm1]-sr)/2.0;
               yi = (hi[enm1][enm1]-si)/2.0;
               csqrt(sqr(yr)-sqr(yi)+xr,2.0*yr*yi+xi,&zzr,&zzi);
               if (yr*zzr+yi*zzi<0.0) {
                  zzr = -zzr;
                  zzi = -zzi;
                  }
               cdiv(xr,xi,yr+zzr,yi+zzi,&xr,&xi);
               sr = sr-xr;
               si = si-xi;
               }
            }
   
         for (i=low;i<=en;i++) {
            hr[i][i] = hr[i][i]-sr;
            hi[i][i] = hi[i][i]-si;
            }
   
         tr = tr+sr;
         ti = ti+si;
         its = its+1;
         itn = itn-1;
         /* reduce to triangle (rows) */
         for (i=l+1;i<=en;i++) {
            sr = hr[i][i-1];
            hr[i][i-1] = 0.0;
            norm = pythag(pythag(hr[i-1][i-1],hi[i-1][i-1]),sr);
            xr = hr[i-1][i-1]/norm;
            wr[i-1] = xr;
            xi = hi[i-1][i-1]/norm;
            wi[i-1] = xi;
            hr[i-1][i-1] = norm;
            hi[i-1][i-1] = 0.0;
            hi[i][i-1] = sr/norm;
   
            for (j=i;j<=n;j++) {
               yr = hr[i-1][j];
               yi = hi[i-1][j];
               zzr = hr[i][j];
               zzi = hi[i][j];
               hr[i-1][j] = xr*yr+xi*yi+hi[i][i-1]*zzr;
               hi[i-1][j] = xr*yi-xi*yr+hi[i][i-1]*zzi;
               hr[i][j] = xr*zzr-xi*zzi-hi[i][i-1]*yr;
               hi[i][j] = xr*zzi+xi*zzr-hi[i][i-1]*yi;
               }
            }
   
         si = hi[en][en];
         if (si!=0.0) {
            norm = pythag(hr[en][en],si);
            sr = hr[en][en]/norm;
            si = si/norm;
            hr[en][en] = norm;
            hi[en][en] = 0.0;
            if (en!=n) {
               for (j=en+1;j<=n;j++) {
                  yr = hr[en][j];
                  yi = hi[en][j];
                  hr[en][j] = sr*yr+si*yi;
                  hi[en][j] = sr*yi-si*yr;
                  }
               }
            }
         /* inverse operation (columns) */
         for (j=l+1;j<=en;j++) {
            xr = wr[j-1];
            xi = wi[j-1];
   
            for (i=1;i<=j;i++) {
               yr = hr[i][j-1];
               yi = 0.0;
               zzr = hr[i][j];
               zzi = hi[i][j];
               if (i!=j) {
                  yi = hi[i][j-1];
                  hi[i][j-1] = xr*yi+xi*yr+hi[j][j-1]*zzi;
                  }
               hr[i][j-1] = xr*yr-xi*yi+hi[j][j-1]*zzr;
               hr[i][j] = xr*zzr+xi*zzi-hi[j][j-1]*yr;
               hi[i][j] = xr*zzi-xi*zzr-hi[j][j-1]*yi;
               }
   
            for (i=low;i<=high;i++) {
               yr = zr[i][j-1];
               yi = zi[i][j-1];
               zzr = zr[i][j];
               zzi = zi[i][j];
               zr[i][j-1] = xr*yr-xi*yi+hi[j][j-1]*zzr;
               zi[i][j-1] = xr*yi+xi*yr+hi[j][j-1]*zzi;
               zr[i][j] = xr*zzr+xi*zzi-hi[j][j-1]*yr;
               zi[i][j] = xr*zzi-xi*zzr-hi[j][j-1]*yi;
               }
            }
   
         if (si!=0.0) {
            for (i=1;i<=en;i++) {
               yr = hr[i][en];
               yi = hi[i][en];
               hr[i][en] = sr*yr-si*yi;
               hi[i][en] = sr*yi+si*yr;
               }
            for (i=low;i<=high;i++) {
               yr = zr[i][en];
               yi = zi[i][en];
               zr[i][en] = sr*yr-si*yi;
               zi[i][en] = sr*yi+si*yr;
               }
            }
         }
      if (itn==0) break;
   
      /* a root found */
      hr[en][en] += tr;
      wr[en] = hr[en][en];
      hi[en][en] += ti;
      wi[en] = hi[en][en];
      en = enm1;
      }
   if (en>=low) /* no convergence after 30*n iterations */
      myerror("COMQR2: no convergence after 30*n iterations");
   
   /* otherwise all roots found, find vectors of upper triangular form */
   norm = 0.0;
   
   for (i=1;i<=n;i++)
      for (j=i;j<=n;j++) {
         tr = dabs(hr[i][j])+dabs(hi[i][j]);
         if (tr>norm)
            norm = tr;
         }
   if (n!=1 && norm!=0.0) {
      for (en=n;en>=2;en--) {
         xr = wr[en];
         xi = wi[en];
         hr[en][en] = 1.0;
         hi[en][en] = 0.0;
         enm1 = en-1;
         for (i=en-1;i>=1;i--) {
            zzr = 0.0;
            zzi = 0.0;
            for (j=i+1;j<=en;j++) {
               zzr = zzr+hr[i][j]*hr[j][en]-hi[i][j]*hi[j][en];
               zzi = zzi+hr[i][j]*hi[j][en]+hi[i][j]*hr[j][en];
               }
   
            yr = xr-wr[i];
            yi = xi-wi[i];
            if (yr==0.0&&yi==0.0) {
               tst1 = norm;
               tst2 = norm+yr;
               yr = tst1;
               while (tst2>tst1) {
                  yr = 0.01*yr;
                  tst2 = norm+yr;
                  }
               }
            cdiv(zzr,zzi,yr,yi,&hr[i][en],&hi[i][en]);
            /* overflow control */
            tr = dabs(hr[i][en])+dabs(hi[i][en]);
            if (tr!=0.0) {
               tst1 = tr;
               tst2 = tst1+1.0/tst1;
               if (tst2<=tst1)
                  for (j=i;j<=en;j++) {
                     hr[j][en] = hr[j][en]/tr;
                     hi[j][en] = hi[j][en]/tr;
                     }
               }
            }
         }
      /* end backsubstitution,  vectors of isolated roots */
      for (i=1;i<=n;i++) 
         if (i<low||i>high) for (j=i;j<=n;j++) { 
               zr[i][j] = hr[i][j];
               zi[i][j] = hi[i][j];
               }
      /* multiply by transform matrix to give vectors of original full matrix */
      for (j=n;j>=low;j--) { 
         m = min(j,high);
   
         for (i=low;i<=high;i++) { 
            zzr = 0.0;
            zzi = 0.0;
   
            for (k=low;k<=m;k++) { 
               zzr = zzr+zr[i][k]*hr[k][j]-zi[i][k]*hi[k][j];
               zzi = zzi+zr[i][k]*hi[k][j]+zi[i][k]*hr[k][j];
               }
   
            zr[i][j] = zzr;
            zi[i][j] = zzi;
            }
         }
      }
}

corth(int n,int low,int high,
      float **ar, float **ai, float *ortr,float *orti)
/*
      given a complex general matrix, this subroutine
      reduces a submatrix situated in rows and columns
      low through high to upper hessenberg form by
      unitary similarity transformations.
 
      on input
         n is the order of the matrix.
         low and high are integers determined by the balancing
           subroutine  cbal.  if  cbal  has not been used,
           set low=1, high=n.
         ar and ai contain the real and imaginary parts,
           respectively, of the complex input matrix.
 
      on output
         ar and ai contain the real and imaginary parts,
           respectively, of the hessenberg matrix.  information
           about the unitary transformations used in the reduction
           is stored in the remaining triangles under the
           hessenberg matrix.
         ortr and orti contain further information about the
           transformations.  only elements low through high are used.
*/
{
   int i,j,m;
   float f,g,h,fi,fr,scale;

   if (high > low+1) for (m=low+1; m<=high-1; m++) {
      h = 0.0;
      ortr[m] = 0.0;
      orti[m] = 0.0;
      scale = 0.0;
      /* scale column (algol tol then not needed) */
      for (i=m; i<=high; i++)
         scale += dabs(ar[i][m-1])+dabs(ai[i][m-1]);

      if (scale!=0.0) {
         for (i=high; i>=m; i--) {
            ortr[i] = ar[i][m-1]/scale;
            orti[i] = ai[i][m-1]/scale;
            h += ortr[i]*ortr[i]+orti[i]*orti[i];
            }

         g = sqrt(h);
         f = pythag(ortr[m],orti[m]);
         if (f==0.0) {
            ortr[m] = g;
            ar[m][m-1] = scale;
         } else {
            h += f*g;
            g /= f;
            ortr[m] *= (1.0+g);
            orti[m] *= (1.0+g);
            }
         /* form (i-(u*ut)/h) * a */
         for (j=m;j<=n;j++) {
            fr = 0.0;
            fi = 0.0;
            /* for i=high step -1 until m do -- */
            for (i=high;i>=m;i--) {
               fr += ortr[i]*ar[i][j]+orti[i]*ai[i][j];
               fi += ortr[i]*ai[i][j]-orti[i]*ar[i][j];
               }

            fr /= h;
            fi /= h;

            for (i=m;i<=high;i++) {
               ar[i][j] += -fr*ortr[i]+fi*orti[i];
               ai[i][j] += -fr*orti[i]-fi*ortr[i];
               }
            }
         /* form (i-(u*ut)/h)*a*(i-(u*ut)/h) */
         for (i=1;i<=high;i++) {
            fr = 0.0;
            fi = 0.0;
            /* for j=high step -1 until m do -- */
            for (j=high;j>=m;j--) {
               fr = fr+ortr[j]*ar[i][j]-orti[j]*ai[i][j];
               fi = fi+ortr[j]*ai[i][j]+orti[j]*ar[i][j];
               }

            fr /= h;
            fi /= h;

            for (j=m;j<=high;j++) {
               ar[i][j] = ar[i][j]-fr*ortr[j]-fi*orti[j];
               ai[i][j] = ai[i][j]+fr*orti[j]-fi*ortr[j];
               }
            }

         ortr[m] *= scale;
         orti[m] *= scale;
         ar[m][m-1] *= -g;
         ai[m][m-1] *= -g;
         }
      }
}

cbabk2(int n,int low,int high, float *scale,int m, float **zr, float **zi)
/*
      this subroutine forms the eigenvectors of a complex general
      matrix by back transforming those of the corresponding
      balanced matrix determined by  cbal.
 
      on input
         n is the order of the matrix.
         low and high are integers determined by  cbal.
         scale contains information determining the permutations
           and scaling factors used by  cbal.
         m is the number of eigenvectors to be back transformed.
         zr and zi contain the real and imaginary parts,
           respectively, of the eigenvectors to be
           back transformed in their first m columns.
 
      on output
         zr and zi contain the real and imaginary parts,
           respectively, of the transformed eigenvectors
           in their first m columns.
*/
{
   int i,j,k,ii;
   float s;
   
   if (m!=0) {
      if (high!=low)
   
         for (i=low;i<=high;i++) {
            s = scale[i];
            /* left hand eigenvectors are back transformed if the 
               foregoing statement is replaced by s=1.0/scale[i] */
            for (j=1;j<=m;j++) {
               zr[i][j] = zr[i][j]*s;
               zi[i][j] = zi[i][j]*s;
               }
            }
      for (ii=1;ii<=n;ii++) {
         i = ii;
         if (i<low || i>high) {
            if (i<low) i = low-ii;
            k = scale[i];
            if (k!=i)
   
               for (j=1;j<=m;j++) {
                  mswap(zr,i,j,k,j);
                  mswap(zi,i,j,k,j);
                  }
            }
         }
      }
}

cbal(int n, float **ar, float **ai,int *low,int *high, float *scale)
/*
      this subroutine balances a complex matrix and isolates
      eigenvalues whenever possible.
 
      on input
         n is the order of the matrix.
         ar and ai contain the real and imaginary parts,
           respectively, of the complex matrix to be balanced.
 
      on output
         ar and ai contain the real and imaginary parts,
           respectively, of the balanced matrix.
         low and high are two integers such that ar[i][j] and ai[i][j]
           are equal to zero if
            (1) i is greater than j and
            (2) j=1,...,low-1 or i=high+1,...,n.
         scale contains information determining the
            permutations and scaling factors used.
 
      suppose that the principal submatrix in rows low through high
      has been balanced, that p[j] denotes the index interchanged
      with j during the permutation step, and that the elements
      of the diagonal matrix used are denoted by d[i][j].  then
         scale[j] = p[j],    for j = 1,...,low-1
                  = d[j][j]      j = low,...,high
                  = p[j]         j = high+1,...,n.
      the order in which the interchanges are made is n to high+1,
      then 1 to low-1.
 
      note that 1 is returned for high if high is zero formally.
 
      the algol procedure exc contained in cbalance appears in
      cbal  in line.  (note that the algol roles of identifiers
      k,l have been reversed.)
*/
{
   int i,j,k,l,m,zflg1,zflg2;
   int noconv;
   float c,f,g,r,s,b2,radix;
   
   radix = 16.0;
   
   b2 = sqr(radix);
   k = 1;
   l = n;
   for (;;) {
      for (j=l;j>=1;j--) {
         zflg1 = 0;
         for (i=1;i<=l;i++) 
            if ( i!=j && (ar[j][i]!=0.0 || ai[j][i]!=0.0) ) zflg1++;
         if (zflg1==0) break;
         }

      if (zflg1==0) {
         m = l;
         /* in-line procedure for row and column exchange */
         scale[m] = j;
         if (j!=m) {
            for (i=1;i<=l;i++) {
               mswap(ar,i,j,i,m);
               mswap(ai,i,j,i,m);
               }
            for (i=k;i<=n;i++) {
               mswap(ar,j,i,m,i);
               mswap(ai,j,i,m,i);
               }
            }
         /* search for rows isolating an eigenvalue and push them down */
         if (l>1) l--;
      } else {
         for (;;) {
            /* search for columns isolating an eigenvalue and push them left */
            for (j=k;j<=l;j++) {
               zflg2=0;
               for (i=k;i<=l;i++) 
                  if ( i!=j && (ar[i][j]!=0.0 || ai[i][j]!=0.0) ) zflg2++;
               if (zflg2==0) break;
               }
            if (zflg2>0) break;
            else {
               m = k;
               /* in-line procedure for row and column exchange */
               scale[m] = j;
               if (j!=m) {
                  for (i=1;i<=l;i++) {
                     mswap(ar,i,j,i,m);
                     mswap(ai,i,j,i,m);
                     }
                  for (i=k;i<=n;i++) {
                     mswap(ar,j,i,m,i);
                     mswap(ai,j,i,m,i);
                     }
                  }
               k = k+1;
               }
            }
         }
      if (zflg1>0 && zflg2>0) break;
      else if (zflg1==0 && l==1) break;
      }

   /* now balance the submatrix in rows k to l */
   if (zflg1>0) {
      for (i=k;i<=l;i++) scale[i] = 1.0;
      noconv = TRUE;
      while (noconv) { /* iterative loop for norm reduction */
         noconv = FALSE;
         for (i=k;i<=l;i++) {
            c = 0.0;
            r = 0.0;
      
            for (j=k;j<=l;j++) if (j!=i) {
               c += dabs(ar[j][i])+dabs(ai[j][i]);
               r += dabs(ar[i][j])+dabs(ai[i][j]);
               }
            if (c!=0.0 && r!=0.0) { /* guard against c,r = 0 from underflow */
               g = r/radix;
               f = 1.0;
               s = c+r;
               while (c<g) {
                  f *= radix;
                  c *= b2;
                  }
               g = r*radix;
               while (c>=g) {
                  f /= radix;
                  c /= b2;
                  }
               if ((c+r)/f<0.95*s) {      /* now balance */
                  g = 1.0/f;
                  scale[i] *= f;
                  noconv = TRUE;
      
                  for (j=k;j<=n;j++) {
                     ar[i][j] *= g;
                     ai[i][j] *= g;
                     }
      
                  for (j=1;j<=l;j++) {
                     ar[j][i] *= f;
                     ai[j][i] *= f;
                     }
                  }
               }
            }
         }
      }
   *low = k;
   *high = l;
}

cg(float **ar,float **ai, int n,
   float *wr, float *wi,float **zr,float **zi)
/*
      on input
         n  is the order of the matrix  a=(ar,ai).
         ar  and  ai  contain the real and imaginary parts,
            respectively, of the complex general matrix.
 
      on output
         wr  and  wi  contain the real and imaginary parts,
            respectively, of the eigenvalues.
         zr  and  zi  contain the real and imaginary parts,
            respectively, of the eigenvectors 

         fv1, fv2, and  fv3  are temporary storage arrays.
*/
{
   int is1,is2;
   float *fv1, *fv2, *fv3;

   fv1 = fvector(1,n);
   fv2 = fvector(1,n);
   fv3 = fvector(1,n);

   cbal(n,ar,ai,&is1,&is2,fv1);
   /*  is1 and is2 are integers determined by the balancing
       subroutine  cbal.  if  cbal  has not been used, set is1=1, is2=n.
   */
   corth(n,is1,is2,ar,ai,fv2,fv3);

   /* find both eigenvalues and eigenvectors */
   comqr2(n,is1,is2,fv2,fv3,ar,ai,wr,wi,zr,zi);
   
   cbabk2(n,is1,is2,fv1,n,zr,zi);

   free_fvector(fv1,1,n);
   free_fvector(fv2,1,n);
   free_fvector(fv3,1,n);
}
