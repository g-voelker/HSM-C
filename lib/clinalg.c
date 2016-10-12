#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"
#include "alloc_space.h"
#include "constants.h"
#include "macros.h"
typedef dcomplex **cmatrix_t, *cvector_t;
typedef zomplex **zmatrix_t, *zvector_t;

extern void nrerror();

#define TINY 1.0e-20

void ludcmp(a, n, indx,d)
   int n;
   int *indx;
   dcomplex *d;
   dcomplex **a;
{
   int i,imax,j,k;
   double big,temp;
   dcomplex *vv;
   dcomplex cOne,dum,sum;

   cOne.r=1.0; cOne.i=0.0;

   vv=dcvector(1,n);
   *d=cOne;
   for (i=1;i<=n;i++) {
      big=0.0;
      for (j=1;j<=n;j++)
         if ((temp=Cabs(a[i][j])) > big) big=temp;
      if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
      vv[i]=Cdiv(cOne,Complex(big,0.0));
   }
   for (j=1;j<=n;j++) {
      for (i=1;i<j;i++) {
         sum=a[i][j];
         for (k=1;k<i;k++) sum = Csub(sum,Cmul(a[i][k],a[k][j]));
         a[i][j]=sum;
      }
      big=0.0;
      for (i=j;i<=n;i++) {
         sum=a[i][j];
         for (k=1;k<j;k++)
            sum = Csub(sum,Cmul(a[i][k],a[k][j]));
         a[i][j]=sum;
         dum = Cmul(vv[i],sum);
         if ( (temp=Cabs(dum)) >= big) {
            big=temp;
            imax=i;
         }
      }
      if (j != imax) {
         for (k=1;k<=n;k++) {
            dum=a[imax][k];
            a[imax][k]=a[j][k];
            a[j][k]=dum;
         }
         *d = RCmul(-1.0,(*d));
         vv[imax]=vv[j];
      }
      indx[j]=imax;
      if (Cabs(a[j][j]) == 0.0) a[j][j]=Complex(TINY,TINY);
      if (j != n) {
         dum=Cdiv(cOne,a[j][j]);
         for (i=j+1;i<=n;i++) a[i][j] = Cmul(a[i][j],dum);
      }
   }
   free_dcvector(vv,1,n);
}

#undef TINY

void lubksb(a, n, indx,b)
   int n;
   int *indx;
   dcomplex b[];
   dcomplex **a;
{
   int i,ii=0,ip,j;
   dcomplex sum;

   for (i=1;i<=n;i++) {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii)
         for (j=ii;j<=i-1;j++) sum = Csub(sum,Cmul(a[i][j],b[j]));
      else if (Cabs(sum)) ii=i;
      b[i]=sum;
   }
   for (i=n;i>=1;i--) {
      sum=b[i];
      for (j=i+1;j<=n;j++) sum = Csub(sum,Cmul(a[i][j],b[j]));
      b[i]=Cdiv(sum,a[i][i]);
   }
}

void printCmat(c, n) 
   int n;
   cmatrix_t c;
{
   int i,j;

   for( i = 1; i <= n; i++ ) {
      for( j = 1; j <= n; j++ ) 
         printf( "(%7.4lf,%7.4lf)", c[i][j].r, c[i][j].i ); 
      putchar( '\n' ); 
   } 
}


void Cmatvectmult(a, b,c, n) 
   int n;
   cvector_t b,c;
   cmatrix_t a;
{
/* c = a * b */

   int i, j;
   dcomplex czero;

   czero = Complex(0.0,0.0);
   for( i = 1; i <= n; i++ ) {
      c[i] = czero;   
      for( j = 1; j <= n; j++ )
         c[i] = Cadd(c[i],Cmul(a[i][j],b[j])); 
   } 
}

void Cmatmult(a,b,c, n) 
   int n;
   cmatrix_t a,b,c;
{
/* c = a * b */

   int i, j, k;
   dcomplex czero;

   czero = Complex(0.0,0.0);

   for( i = 1; i <= n; i++ )
      for( k = 1; k <= n; k++ ) {
         c[i][k] = czero;   
         for( j = 1; j <= n; j++ )
            c[i][k] = Cadd(c[i][k],Cmul(a[i][j],b[j][k])); 
      } 
}


void Cinvert(a,b, n) 
   int n;
   cmatrix_t a,b;
{
/* b = 1/a where a is nxn */

   int *indx, i, j;
   dcomplex d, *col, czero, cone;

   czero = Complex(0.0,0.0);
   cone = Complex(1.0,0.0);

   indx = ivector( 1, n );
   col = dcvector( 1, n );
   ludcmp( a, n, indx, &d );
   for( j = 1; j <= n; j++ ) {
      for( i = 1; i <= n; col[i++] = czero );
      col[j] = cone;
      lubksb( a, n, indx, col );
      for( i = 1; i <= n; i++ ) b[i][j] = col[i]; 
   } 
   free_ivector( indx, 1, n );
   free_dcvector( col, 1, n );
}

dcomplex Cdotprod(a,b, n) 
   int n;
   cvector_t a,b;
{
/* a* dot b */

   int i;
   dcomplex dp;

   dp = Cmul(a[1],b[1]);
   for( i = 2; i <= n; i++ ) dp = Cadd(dp,Cmul(a[i],b[i]));
   return(dp);
}

/* ------------------------------------------------------------------
   ZOMPLEX ROUTINES
   ------------------------------------------------------------------
*/

#define TINY 1.0e-20

void zludcmp(a, n, indx,d)
   int n;
   int *indx;
   zomplex *d;
   zomplex **a;
{
   int i,imax,j,k;
   double big,temp;
   zomplex *vv;
   zomplex cOne,dum,sum;

   cOne.re=1.0; cOne.im=0.0;

   vv=zvector(1,n);
   *d=cOne;
   for (i=1;i<=n;i++) {
      big=0.0;
      for (j=1;j<=n;j++)
         if ((temp=Zabs(a[i][j])) > big) big=temp;
      if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
      vv[i]=Zdiv(cOne,Zomplex(big,0.0));
   }
   for (j=1;j<=n;j++) {
      for (i=1;i<j;i++) {
         sum=a[i][j];
         for (k=1;k<i;k++) sum = Zsub(sum,Zmul(a[i][k],a[k][j]));
         a[i][j]=sum;
      }
      big=0.0;
      for (i=j;i<=n;i++) {
         sum=a[i][j];
         for (k=1;k<j;k++)
            sum = Zsub(sum,Zmul(a[i][k],a[k][j]));
         a[i][j]=sum;
         dum = Zmul(vv[i],sum);
         if ( (temp=Zabs(dum)) >= big) {
            big=temp;
            imax=i;
         }
      }
      if (j != imax) {
         for (k=1;k<=n;k++) {
            dum=a[imax][k];
            a[imax][k]=a[j][k];
            a[j][k]=dum;
         }
         *d = RZmul(-1.0,(*d));
         vv[imax]=vv[j];
      }
      indx[j]=imax;
      if (Zabs(a[j][j]) == 0.0) a[j][j]=Zomplex(TINY,TINY);
      if (j != n) {
         dum=Zdiv(cOne,a[j][j]);
         for (i=j+1;i<=n;i++) a[i][j] = Zmul(a[i][j],dum);
      }
   }
   free_zvector(vv,1,n);
}

#undef TINY

void zlubksb(a, n, indx,b)
   int n;
   int *indx;
   zomplex b[];
   zomplex **a;
{
   int i,ii=0,ip,j;
   zomplex sum;

   for (i=1;i<=n;i++) {
      ip=indx[i];
      sum=b[ip];
      b[ip]=b[i];
      if (ii)
         for (j=ii;j<=i-1;j++) sum = Zsub(sum,Zmul(a[i][j],b[j]));
      else if (Zabs(sum)) ii=i;
      b[i]=sum;
   }
   for (i=n;i>=1;i--) {
      sum=b[i];
      for (j=i+1;j<=n;j++) sum = Zsub(sum,Zmul(a[i][j],b[j]));
      b[i]=Zdiv(sum,a[i][i]);
   }
}


void printZmat(c, n) 
   int n;
   zmatrix_t c;
{
   int i,j;

   for( i = 1; i <= n; i++ ) {
      for( j = 1; j <= n; j++ ) 
         printf( "(%7.4lf,%7.4lf)", c[i][j].re, c[i][j].im ); 
      putchar( '\n' ); 
   } 
}


void Zmatvectmult(a, b,c, n) 
   int n;
   zvector_t b,c;
   zmatrix_t a;
{
/* c = a * b */

   int i, j;
   zomplex czero;

   czero = Zomplex(0.0,0.0);
   for( i = 1; i <= n; i++ ) {
      c[i] = czero;   
      for( j = 1; j <= n; j++ )
         c[i] = Zadd(c[i],Zmul(a[i][j],b[j])); 
   } 
}

void Zmatmult(a,b,c, n) 
   int n;
   zmatrix_t a,b,c;
{
/* c = a * b */

   int i, j, k;
   zomplex czero;

   czero = Zomplex(0.0,0.0);

   for( i = 1; i <= n; i++ )
      for( k = 1; k <= n; k++ ) {
         c[i][k] = czero;   
         for( j = 1; j <= n; j++ )
            c[i][k] = Zadd(c[i][k],Zmul(a[i][j],b[j][k])); 
      } 
}


void Zinvert(a,b, n) 
   int n;
   zmatrix_t a,b;
{
/* b = 1/a where a is nxn */

   int *indx, i, j;
   zomplex d, *col, czero, cone;

   czero = Zomplex(0.0,0.0);
   cone = Zomplex(1.0,0.0);

   indx = ivector( 1, n );
   col = zvector( 1, n );
   zludcmp( a, n, indx, &d );
   for( j = 1; j <= n; j++ ) {
      for( i = 1; i <= n; col[i++] = czero );
      col[j] = cone;
      zlubksb( a, n, indx, col );
      for( i = 1; i <= n; i++ ) b[i][j] = col[i]; 
   } 
   free_ivector( indx, 1, n );
   free_zvector( col, 1, n );
}

zomplex Zdotprod(a,b, n) 
   int n;
   zvector_t a,b;
{
/* a* dot b */

   int i;
   zomplex dp;

   dp = Zmul(a[1],b[1]);
   for( i = 2; i <= n; i++ ) dp = Zadd(dp,Zmul(a[i],b[i]));
   return(dp);
}

