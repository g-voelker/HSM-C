#include <stdio.h>
#include <math.h>
#include "complex.h"
#include "alloc_space.h"
#include "macros.h"
typedef double **matrix_t, *vector_t;

#define TINY 1.0e-20;

void ludcmp(a,n,indx,d)
   int n,*indx;
   double **a,*d;
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv,*dvector();
	void nrerror(),free_dvector();

	vv=dvector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=dabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine LUDCMP");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*dabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,1,n);
}

#undef TINY

void lubksb(a,n,indx,b)
   double **a,b[];
   int n,*indx;
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

void printmat(c,n) 
   int n;
   matrix_t c;
{
   int i,j;

   for( i = 1; i <= n; i++ ) {
      for( j = 1; j <= n; j++ ) printf( "%7.4lf", c[i][j] ); 
		putchar( '\n' ); 
   } 
}


void matvectmult(a,b,c, n) 
   int n;
   vector_t b,c;
   matrix_t a;
{
/* c = a * b */

	int i, j;

	for( i = 1; i <= n; i++ ) {
		c[i] = 0.0;	
		for( j = 1; j <= n; j++ )
			c[i] += a[i][j] * b[j]; 
	} 
}

void matmult(a,b,c, n) 
   int n;
   matrix_t a,b,c;
{
/* c = a * b */

	int i, j, k;

	for( i = 1; i <= n; i++ )
		for( k = 1; k <= n; k++ ) {
			c[i][k] = 0.0;	
			for( j = 1; j <= n; j++ )
				c[i][k] += a[i][j] * b[j][k]; 
		} 
}


void invert(a,b, n) 
   int n;
   matrix_t a,b;
{
/* b = 1/a where a is nxn */

   int *indx, i, j;
   double d, *col;

   indx = ivector( 1, n );
   col = dvector( 1, n );
   ludcmp( a, n, indx, &d );
   for( j = 1; j <= n; j++ ) {
      for( i = 1; i <= n; col[i++] = 0.0 );
      col[j] = 1.0;
      lubksb( a, n, indx, col );
      for( i = 1; i <= n; i++ ) b[i][j] = col[i]; 
   } 
   free_ivector( indx, 1, n );
   free_dvector( col, 1, n );
}

