/* ========================================================================
   This file contains various programs to take derivatives of functions of
   one or more variable.  1-dim derivatives suppose points are defined in 
   the index range 1..n.  2 and 3-dim derivatives suppose points are defined 
   in the index range 0..n, 0..m, etc.
 
      intplderiv1 - take first order derivative given two points
      intplderiv - take derivative at interpolated point
      diff       - find derivative of f(x)  using O(x^2) accurate method

      pdiff2DO1   - find partial derivative of f(x,y) wrt x or y
      pdiff2D     - find partial derivative of f(x,y) wrt x or y
      pdiff2D_gen - find partial derivative of f(x,y) wrt variable x or y

    Z/Cpdiff2D   - find partial derivative of complex f(x,y) wrt x or y
    Z/Cpdiff2DO1 - find partial derivative of complex f(x,y) wrt x or y
                   This method is only O(x,y,z) accurate
    Z/Cpdiff2D_gen - find partial derivative of complex f(x,y) wrt var. x or y
    Z/CFpdiff2D      - find partial derivative of complex fourier f(x,y) 
    Z/CFpdiff2D_gen  - find partial derivative of complex fourier f(x,y) 

    Z/Cpdiff3D   - find partial derivative of complex f(x,y,z) wrt x,y,or z.
    Z/Cpdiff3DO1 - find partial derivative of complex f(x,y,z) wrt x,y,or z.
                   This method is only O(x,y,z) accurate
    Z/Cpdiff3D_gen - find partial derivative of complex fourier f(x,y,z) 
                     wrt variable x, y, or z
    Z/CFpdiff3D      - find partial derivative of complex fourier f(x,y,z) 
    Z/CFpdiff3D_gen  - find partial derivative of complex fourier f(x,y,z) 

   ========================================================================
*/
#include <stdio.h>
#include "myerror.h"
#include "complex.h"
#include "alloc_space.h"
#include "constants.h"
#include "macros.h"

/* ==========================================================
   Find O(x^1) derivative 
   ==========================================================
 */
double intplderiv1(x1,x2, y1,y2)
   double x1,x2,y1,y2;
{
  return( (y2-y1)/(x2-x1) ); 
}

/* ==========================================================
   Find derivative at an interpolated point (O(x^2) accurate)
   ==========================================================
 */
double intplderiv(xval, x1,x2,x3, y1,y2,y3)
   double xval, x1,x2,x3,y1,y2,y3;
{
  double xv2,d12,d23,d31,d1,d2,d3;

  xv2=2*xval;

  d12 = x1-x2;  d23 = x2-x3;  d31 = x3 - x1;
  d1 = -y1/(d12*d31);  d2 = -y2/(d23*d12);  d3 = -y3/(d31*d23);

  return( d1*(xv2-x2-x3) + d2*(xv2-x1-x3) + d3*(xv2-x1-x2) ); 
}

/* =========================================================
   Find derivative of a function fn[i](x[i])
   =========================================================
*/
void diff(fn,x,nind, fnp)
  int nind;
  double *fn,*x,*fnp;
{
  int i;

  fnp[1] = intplderiv(x[1], x[1],x[2],x[3],fn[1],fn[2],fn[3]);
  fnp[nind] = intplderiv(x[nind], x[nind],x[nind-1],x[nind-2],
                                  fn[nind],fn[nind-1],fn[nind-2]);
  for (i=2;i<nind;i++) 
     fnp[i] = intplderiv(x[i], x[i-1],x[i],x[i+1],fn[i-1],fn[i],fn[i+1]);
}

/* =========================================================
   DIFFO1: Jan 10/94

   Find derivative of a function fn[i](x[i]) to first order
   =========================================================
*/
void diffO1(fn,x,nind, fnp,isign)
   int nind, isign;
   double *fn,*x,*fnp;
{
  int i;

  if (isign == 1) {
     for (i=2;i<=nind;i++) 
        fnp[i] = (fn[i]-fn[i-1])/(x[i]-x[i-1]);
     fnp[1] = fnp[2];
  } else {
     for (i=1;i<nind;i++) 
        fnp[i] = (fn[i+1]-fn[i])/(x[i+1]-x[i]);
     fnp[nind] = fnp[nind-1];
  }
}

/* =========================================================
   diff1D Find order 2 accurate derivative of a function 
   =========================================================
*/
void diff1D(f, ox,nx, dx, df)
  int ox,nx;
  double dx;
  double *f,*df;
{
  int i;
  double idx2;

  idx2 = 0.5/dx;
  for (i=ox+1;i<nx;i++) df[i] = idx2*(f[i+1]-f[i-1]);
  df[ox] = idx2*(-3*f[ox] + 4*f[ox+1] - f[ox+2]);
  df[nx] = idx2*(3*f[nx] - 4*f[nx-1] + f[nx-2]);
}

/* =========================================================
   Pdiff1D: Find derivative of a periodic function 
   =========================================================
*/
void Pdiff1D(f, ox,nx, dx, df)
  int ox,nx;
  double dx;
  double *f,*df;
{
  int i;
  double idx2;

  idx2 = 0.5/dx;
  for (i=ox+1;i<nx;i++) df[i] = idx2*(f[i+1]-f[i-1]);
  df[ox]=df[nx] = idx2*(f[ox+1]-f[nx-1]);
}

/* =========================================================================
   PDIFF2DO1: Sept 23/93 

   Routines to take find derivative at an interpolated point
   and to take derivative of a function fn[i](x[i]).  Unlike
   diff, pdiff2D and other 2D functions assume that the matrix
   and index functions are defined on a range 0..n (not 1..n)

   Function is assumed to be defined at evenly spaced grid points
   Method is only Order(1) accurate.  The derivative is taken
   with respect to "|index|".  
   If index is positive then (i.e index==1):  df[i][j] = f[i+1][j]-f[i][j]  
   If index is negative then (i.e index==-1): df[i][j] = f[i][j]-f[i-1][j]  

   Note that df is not defined over the full index range of f.  It is up
   to the calling program to take care of boundaries but, by default, 
   df inserts zeros where it is not defined.
   =========================================================================
*/

void pdiff2DO1(f, ox,nx,oy,ny, dxy,index, df)
   int ox,nx,oy,ny,index;
   double dxy;
   double **f,**df;
{
  int aind;
  int i,j,ip1,jp1,im1,jm1;
  double idx,idy;

  aind = iabs(index);
  if (aind==0 || aind>2) myerror("PDIFF2DO1: index out of bounds");

  if (index==1 && ox!=1) myerror("PDIFF2DO1: ox should equal one");
  if (index==2 && oy!=1) myerror("PDIFF2DO1: oy should equal one");
  if (index==-1 && ox!=0) myerror("PDIFF2DO1: ox should equal zero");
  if (index==-2 && oy!=0) myerror("PDIFF2DO1: oy should equal zero");

  if (index==1) {  /* take derivative wrt x, assume evenly spaced */
     idx = 1.0/dxy;
     for (i=1;i<nx;i++) { 
        ip1=i+1;
        for (j=oy;j<=ny;j++) df[i][j] = idx*(f[ip1][j]-f[i][j]);
     }
     for (j=oy;j<=ny;j++) {
        df[0][j] = df[1][j];
        df[nx][j] = df[nx-1][j];
     }
  } else if (index==2) {    /* take derivative wrt y, assume evenly spaced */
     idy = 1.0/dxy;
     for (j=1;j<ny;j++) {
        jp1=j+1;
        for (i=ox;i<=nx;i++) df[i][j] = idy*(f[i][jp1]-f[i][j]);
     }
     for (i=ox;i<=nx;i++) {
        df[i][0] = df[i][1];
        df[i][ny] = df[i][ny-1];
     }
  } else if (index==-1) {  /* take derivative wrt x, assume evenly spaced */
     idx = 1.0/dxy;
     for (i=1;i<=nx;i++) { 
        im1=i-1;
        for (j=oy;j<=ny;j++) df[i][j] = idx*(f[i][j]-f[im1][j]);
     }
  } else if (index==-2) {    /* take derivative wrt y, assume evenly spaced */
     idy = 1.0/dxy;
     for (j=1;j<=ny;j++) {
        jm1=j-1;
        for (i=ox;i<=nx;i++) df[i][j] = idy*(f[i][j]-f[i][jm1]);
     }
  }
}

/* =========================================================================
   PPDIFF2DO1: Dec 6/93

   Routines to take find derivative at an interpolated point
   and to take derivative of a function fn[i](x[i]).  This function is
   similar to pdiff2DO1 except for the treatment of boundaries which
   are assumed to be periodic.

   Function is assumed to be defined at evenly spaced grid points
   Method is only Order(1) accurate.  The derivative is taken
   with respect to "|index|".  
   If index is positive then (i.e index==1):  df[i][j] = f[i+1][j]-f[i][j]  
   If index is negative then (i.e index==-1): df[i][j] = f[i][j]-f[i-1][j]  

   Note that df is not defined over the full index range of f.  It is up
   to the calling program to take care of boundaries but, by default, 
   df inserts zeros where it is not defined.
   =========================================================================
*/

void Ppdiff2DO1(f, ox,nx,oy,ny, dxy,index, df)
   int ox,nx,oy,ny,index;
   double dxy;
   double **f,**df;
{
  int aind;
  int i,j,ip1,jp1,im1,jm1;
  double idx,idy;

  aind = iabs(index);
  if (aind==0 || aind>2) myerror("PPDIFF2DO1: index out of bounds");

  if (aind==1 && ox!=0) myerror("PPDIFF2DO1: ox should equal zero");
  if (aind==2 && oy!=0) myerror("PPDIFF2DO1: oy should equal zero");

  if (index==1) {  /* take derivative wrt x, assume evenly spaced */
     idx = 1.0/dxy;
     for (i=0;i<nx;i++) { 
        ip1=i+1;
        for (j=oy;j<=ny;j++) df[i][j] = idx*(f[ip1][j]-f[i][j]);
     }
     for (j=oy;j<=ny;j++) df[nx][j] = df[0][j];
  } else if (index==2) {    /* take derivative wrt y, assume evenly spaced */
     idy = 1.0/dxy;
     for (j=0;j<ny;j++) {
        jp1=j+1;
        for (i=ox;i<=nx;i++) df[i][j] = idy*(f[i][jp1]-f[i][j]);
     }
     for (i=ox;i<=nx;i++) df[i][ny] = df[i][0];
  } else if (index==-1) {  /* take derivative wrt x, assume evenly spaced */
     idx = 1.0/dxy;
     for (i=1;i<=nx;i++) { 
        im1=i-1;
        for (j=oy;j<=ny;j++) df[i][j] = idx*(f[i][j]-f[im1][j]);
     }
     for (j=oy;j<=ny;j++) df[0][j] = df[nx][j];
  } else if (index==-2) {    /* take derivative wrt y, assume evenly spaced */
     idy = 1.0/dxy;
     for (j=1;j<=ny;j++) {
        jm1=j-1;
        for (i=ox;i<=nx;i++) df[i][j] = idy*(f[i][j]-f[i][jm1]);
     }
     for (i=ox;i<=nx;i++) df[i][0] = df[i][ny];
  }
}

/* =========================================================================
   ZPDIFF2DO1: Feb 16/94

   Routines to take find derivative at an interpolated point
   and to take derivative of a function fn[i](x[i]).  This function is
   similar to pdiff2DO1 except for the treatment of boundaries which
   are assumed to be zero.

   Function is assumed to be defined at evenly spaced grid points
   Method is only Order(1) accurate.  The derivative is taken
   with respect to "|index|".  
   If index is positive then (i.e index==1):  df[i][j] = f[i+1][j]-f[i][j]  
   If index is negative then (i.e index==-1): df[i][j] = f[i][j]-f[i-1][j]  

   Note that df is not defined over the full index range of f.  It is up
   to the calling program to take care of boundaries but, by default, 
   df inserts zeros where it is not defined.
   =========================================================================
*/

void Zpdiff2DO1(f, ox,nx,oy,ny, dxy,index, df)
   int ox,nx,oy,ny,index;
   double dxy;
   double **f,**df;
{
  int aind;
  int i,j,ip1,jp1,im1,jm1;
  double idx,idy;

  aind = iabs(index);
  if (aind==0 || aind>2) myerror("ZPDIFF2DO1: index out of bounds");

  if (index==1 && ox!=1) myerror("ZPDIFF2DO1: ox should equal one");
  if (index==2 && oy!=1) myerror("ZPDIFF2DO1: oy should equal one");
  if (index==-1 && ox!=0) myerror("ZPDIFF2DO1: ox should equal zero");
  if (index==-2 && oy!=0) myerror("ZPDIFF2DO1: oy should equal zero");

  if (index==1) {  /* take derivative wrt x, assume evenly spaced */
     idx = 1.0/dxy;
     for (i=1;i<nx;i++) { 
        ip1=i+1;
        for (j=oy;j<=ny;j++) df[i][j] = idx*(f[ip1][j]-f[i][j]);
     }
     for (j=oy;j<=ny;j++) df[nx][j] = df[0][j] = 0.0;
  } else if (index==2) {    /* take derivative wrt y, assume evenly spaced */
     idy = 1.0/dxy;
     for (j=1;j<ny;j++) {
        jp1=j+1;
        for (i=ox;i<=nx;i++) df[i][j] = idy*(f[i][jp1]-f[i][j]);
     }
     for (i=ox;i<=nx;i++) df[i][ny] = df[i][0] = 0.0;
  } else if (index==-1) {  /* take derivative wrt x, assume evenly spaced */
     idx = 1.0/dxy;
     for (i=1;i<=nx;i++) { 
        im1=i-1;
        for (j=oy;j<=ny;j++) df[i][j] = idx*(f[i][j]-f[im1][j]);
     }
  } else if (index==-2) {    /* take derivative wrt y, assume evenly spaced */
     idy = 1.0/dxy;
     for (j=1;j<=ny;j++) {
        jm1=j-1;
        for (i=ox;i<=nx;i++) df[i][j] = idy*(f[i][j]-f[i][jm1]);
     }
  }
}

/* =======================================================================
   PDIFF2D: Sept 23/93 

   Routines to take find derivative at an interpolated point
   and to take derivative of a function fn[i](x[i]).   Method is second order 
   accurate. Unlike diff, pdiff2D and other 2D functions assume that the 
   matrix and index functions are defined on a range 0..n (not 1..n)

   Function is assumed to be defined at evenly spaced grid points
   =======================================================================
*/
void pdiff2D(f, ox,nx,oy,ny, dxy,index, df)
   int ox,nx,oy,ny,index;
   double dxy;
   double **f,**df;
{
  int i,j,im1,ip1,jm1,jp1,Nm1,Nm2;
  double idx2,idy2;

  if (index<1 || index>2) {
     printf("ox=%d, nx=%d,  oy=%d,ny=%d,  dxy=%10.3e, index=%d",
        ox,nx,oy,ny,dxy,index);
     fflush(NULL);
     myerror("PDIFF2D: index out of bounds");
  }

  if (index==1) {  /* take derivative wrt x, assume evenly spaced */
     idx2 = 0.5/dxy;
     for (i=ox+1;i<nx;i++) { 
        im1=i-1;  ip1=i+1;
        for (j=oy;j<=ny;j++) df[i][j] = idx2*(f[ip1][j]-f[im1][j]);
     }
     Nm1=nx-1; Nm2=nx-2;
     for (j=oy;j<=ny;j++) {
        df[ox][j] = idx2*(-3*f[ox][j] + 4*f[ox+1][j] - f[ox+2][j]);
        df[nx][j] = idx2*(3*f[nx][j] - 4*f[Nm1][j] + f[Nm2][j]);
     }
  } else if (index==2) {    /* take derivative wrt y, assume evenly spaced */
     idy2 = 0.5/dxy;
     Nm1=ny-1;  Nm2=ny-2;
     for (i=ox;i<=nx;i++) {
        df[i][oy] = idy2*(-3*f[i][oy] + 4*f[i][oy+1] - f[i][oy+2]);
        df[i][ny] = idy2*(3*f[i][ny] - 4*f[i][Nm1] + f[i][Nm2]);
     }
     for (j=oy+1;j<ny;j++) {
        jm1=j-1;  jp1=j+1;
        for (i=ox;i<=nx;i++) df[i][j] = idy2*(f[i][jp1]-f[i][jm1]);
     }
  }
}

/* =========================================================
   PDIFF2D_gen: Sept 23/93 

   Routines to take find derivative at an interpolated point
   and to take derivative of a function fn[i](x[i]).  Unlike
   diff, pdiff2D and other 2D functions assume that the matrix
   and index functions are defined on a range 0..n (not 1..n)

   Function may be defined at unevenly spaced grid points
   =========================================================
*/

void pdiff2D_gen(f, x,ox,nx, y,oy,ny, index, df)
   int ox,nx,oy,ny,index;
   double *x,*y;
   double **f,**df;
{
  int i,j,im1,ip1,jm1,jp1,Nm1,Nm2;

  if (index<1 || index>2) myerror("PDIFF2D_gen: index out of bounds");

  if (index==1) {    /* take derivative wrt x */
     Nm1=nx-1;  Nm2=nx-2;
     for (j=oy;j<=ny;j++) {
        df[ox][j] = intplderiv(x[ox], x[ox],x[ox+1],x[ox+2],
           f[ox][j],f[ox+1][j],f[ox+2][j]);   
        df[nx][j] = intplderiv(x[nx], x[nx],x[Nm1],x[Nm2],
           f[nx][j],f[Nm1][j],f[Nm2][j]);   
     }
     for (i=ox+1;i<nx;i++) { 
        im1=i-1;  ip1=i+1;
        for (j=oy;j<=ny;j++)
           df[i][j] = intplderiv(x[i], x[im1],x[i],x[ip1],
                                       f[im1][j],f[i][j],f[ip1][j]);
     }
  } else if (index==2) {  /* take derivative wrt y */
     Nm1=ny-1;  Nm2=ny-2;
     for (i=ox;i<=nx;i++) {
        df[i][oy] = intplderiv(y[oy], y[oy],y[oy+1],y[oy+2],
           f[i][oy],f[i][oy+1],f[i][oy+2]);   
        df[i][ny] = intplderiv(y[ny], y[ny],y[Nm1],y[Nm2],
           f[i][ny],f[i][Nm1],f[i][Nm2]);   
     }
     for (j=oy+1;j<ny;j++) {
        jm1=j-1;  jp1=j+1;
        for (i=ox;i<=nx;i++)
           df[i][j] = intplderiv(y[j], y[jm1],y[j],y[jp1],
                                       f[i][jm1],f[i][j],f[i][jp1]);
     }
  }
}

/* =======================================================================
   PPDIFF2D: May 1/94

   Routines to take find derivative at an interpolated point
   and to take derivative of a function fn[i](x[i]).   Method is second order 
   accurate and assumes periodicity in direction of differentiation.

   Function is assumed to be defined at evenly spaced grid points
   =======================================================================
*/
void Ppdiff2D(f,ox,nx,oy,ny, dxy,index, df)
   int ox,nx,oy,ny,index;
   double dxy;
   double **f,**df;
{
  int i,j,im1,ip1,jm1,jp1;
  double idx2,idy2;

  if (index<=0 || index>2) myerror("PPDIFF2D: index out of bounds");

  if (index==1 && ox!=0) myerror("PPDIFF2D: ox should equal zero");
  if (index==2 && oy!=0) myerror("PPDIFF2D: oy should equal zero");

  if (index==1) {  /* take derivative wrt x, assume evenly spaced */
     idx2 = 0.5/dxy;
     for (i=ox+1;i<nx;i++) { 
        im1=i-1;  ip1=i+1;
        for (j=oy;j<=ny;j++) df[i][j] = idx2*(f[ip1][j]-f[im1][j]);
     }
     for (j=oy;j<=ny;j++) 
        df[nx][j]=df[ox][j] = idx2*(f[nx-1][j]-f[ox+1][j]);

  } else if (index==2) {    /* take derivative wrt y, assume evenly spaced */
     idy2 = 0.5/dxy;
     for (j=oy+1;j<ny;j++) {
        jm1=j-1;  jp1=j+1;
        for (i=ox;i<=nx;i++) df[i][j] = idy2*(f[i][jp1]-f[i][jm1]);
     }
     for (i=ox;i<=nx;i++) 
        df[i][ny]=df[i][oy] = idy2*(f[i][ny-1]-f[i][oy+1]);
  }
}

/* ============================================================================
   CPDIFF2DO1:  Sept 20/93

   Find partial derivative of a complex 2D matrix with respect to index.
   using a _first_ order finite difference scheme.  Derivative is presumed
   to be calculated at the intermediate gridpoint.  Direction of the
   index depends on the sign of index (i.e. see comments in PDIFF2DO1).
   By default, where df is not defined zeros are inserted.
   ============================================================================
*/
void Cpdiff2DO1(f,ox,nx,oy,ny, dxy,index, df)
   int ox,nx,oy,ny,index;
   double dxy;
   dcomplex **f,**df;
{
  int aind;
  int i,j,ip1,jp1,im1,jm1;
  double idx,idy;

  aind = iabs(index);
  if (aind==0 || aind>2) myerror("CPDIFF2DO1: index out of bounds");

  if (index==1 && ox!=1) myerror("CPDIFF2DO1: ox should equal one");
  if (index==2 && oy!=1) myerror("CPDIFF2DO1: oy should equal one");
  if (index==-1 && ox!=0) myerror("CPDIFF2DO1: ox should equal zero");
  if (index==-2 && oy!=0) myerror("CPDIFF2DO1: oy should equal zero");

  if (index==1) {
     idx = 1.0/dxy;
     for (i=1;i<nx;i++) {
        ip1=i+1;
        for (j=oy;j<=ny;j++) {
           df[i][j].r = idx*(f[ip1][j].r-f[i][j].r);
           df[i][j].i = idx*(f[ip1][j].i-f[i][j].i);
        }
     }
     for (j=oy;j<=ny;j++) {
        df[0][j].r = df[1][j].r;
        df[0][j].i = df[1][j].i;
        df[nx][j].r = df[nx-1][j].r;
        df[nx][j].i = df[nx-1][j].i;
     }
  } else if (index==2) {
     idy = 1.0/dxy;
     for (j=1;j<ny;j++) { 
        jp1=j+1;
        for (i=ox;i<=nx;i++)  { 
           df[i][j].r = idy*(f[i][jp1].r-f[i][j].r); 
           df[i][j].i = idy*(f[i][jp1].i-f[i][j].i);
        }
     }
     for (i=ox;i<=nx;i++) {
        df[i][0].r = df[i][1].r;
        df[i][0].i = df[i][1].i;
        df[i][ny].r = df[i][ny-1].r;
        df[i][ny].i = df[i][ny-1].i;
     }
  } else if (index==-1) {
     idx = 1.0/dxy;
     for (i=1;i<=nx;i++) {
        im1=i-1;
        for (j=oy;j<=ny;j++) {
           df[i][j].r = idx*(f[i][j].r-f[im1][j].r);
           df[i][j].i = idx*(f[i][j].i-f[im1][j].i);
        }
     }
  } else if (index==-2) {
     idy = 1.0/dxy;
     for (j=1;j<=ny;j++) { 
        jm1=j-1;
        for (i=ox;i<=nx;i++)  { 
           df[i][j].r = idy*(f[i][j].r-f[i][jm1].r); 
           df[i][j].i = idy*(f[i][j].i-f[i][jm1].i);
        }
     }
  }
}

/* ============================================================================
   CPPDIFF2DO1:  Dec 7/93

   Find partial derivative of a complex 2D matrix with respect to index.
   using a _first_ order finite difference scheme.  Derivative is presumed
   to be calculated at the intermediate gridpoint.  Direction of the
   index depends on the sign of index (i.e. see comments in PDIFF2DO1).
   By default, where df is not defined zeros are inserted.

   This procedure differs from CPDIFF2DO1 in that the domain is assumed to
   be periodic at the boundary.
   ============================================================================
*/
void CPpdiff2DO1(f,ox,nx,oy,ny, dxy,index, df)
   int ox,nx,oy,ny,index;
   double dxy;
   dcomplex **f,**df;
{
  int aind;
  int i,j,ip1,jp1,im1,jm1;
  double idx,idy;

  aind = iabs(index);
  if (aind==0 || aind>2) myerror("CPPDIFF2DO1: index out of bounds");

  if (aind==1 && ox!=0) myerror("CPPDIFF2DO1: ox should equal zero");
  if (aind==2 && oy!=0) myerror("CPPDIFF2DO1: oy should equal zero");

  if (index==1) {
     idx = 1.0/dxy;
     for (i=0;i<nx;i++) {
        ip1=i+1;
        for (j=oy;j<=ny;j++) {
           df[i][j].r = idx*(f[ip1][j].r-f[i][j].r);
           df[i][j].i = idx*(f[ip1][j].i-f[i][j].i);
        }
     }
     for (j=oy;j<=ny;j++) df[nx][j] = df[0][j];
  } else if (index==2) {
     idy = 1.0/dxy;
     for (j=0;j<ny;j++) { 
        jp1=j+1;
        for (i=ox;i<=nx;i++)  { 
           df[i][j].r = idy*(f[i][jp1].r-f[i][j].r); 
           df[i][j].i = idy*(f[i][jp1].i-f[i][j].i);
        }
     }
     for (i=ox;i<=nx;i++) df[i][ny] = df[i][0];
  } else if (index==-1) {
     idx = 1.0/dxy;
     for (i=1;i<=nx;i++) {
        im1=i-1;
        for (j=oy;j<=ny;j++) {
           df[i][j].r = idx*(f[i][j].r-f[im1][j].r);
           df[i][j].i = idx*(f[i][j].i-f[im1][j].i);
        }
     }
     for (j=oy;j<=ny;j++) df[0][j] = df[nx][j];
  } else if (index==-2) {
     idy = 1.0/dxy;
     for (j=1;j<=ny;j++) { 
        jm1=j-1;
        for (i=ox;i<=nx;i++)  { 
           df[i][j].r = idy*(f[i][j].r-f[i][jm1].r); 
           df[i][j].i = idy*(f[i][j].i-f[i][jm1].i);
        }
     }
     for (i=ox;i<=nx;i++) df[i][0] = df[i][ny];
  }
}

/* ============================================================================
   CPDIFF2D:  Sept 5/93

   Find partial derivative of a complex 2D matrix with respect to index.
   using a second order finite difference scheme.  Function is defined
   at evenly spaced grid points with spacing dxy in given direction of
   differentiation.
   ============================================================================
*/
void Cpdiff2D(f,ox,nx,oy,ny, dxy,index, df)
   int ox,nx,oy,ny,index;
   double dxy;
   dcomplex **f,**df;
{
  int i,j,Nm1,Nm2,im1,ip1,jm1,jp1;
  double idx2,idy2;

  if (index<1 || index>2) myerror("CPDIFF2D: index out of bounds");

  if (index==1) {
     idx2 = 0.5/dxy;
     for (i=ox+1;i<nx;i++) {
        im1=i-1;  ip1=i+1;
        for (j=oy;j<=ny;j++) {
           df[i][j].r = idx2*(f[ip1][j].r-f[im1][j].r);
           df[i][j].i = idx2*(f[ip1][j].i-f[im1][j].i);
        }
     }
     Nm1=nx-1;  Nm2=nx-2;
     for (j=oy;j<=ny;j++) {
        df[ox][j].r = idx2*(-3*f[ox][j].r + 4*f[ox+1][j].r - f[ox+2][j].r);
        df[ox][j].i = idx2*(-3*f[ox][j].i + 4*f[ox+1][j].i - f[ox+2][j].i);
        df[nx][j].r = idx2*(3*f[nx][j].r - 4*f[Nm1][j].r + f[Nm2][j].r);
        df[nx][j].i = idx2*(3*f[nx][j].i - 4*f[Nm1][j].i + f[Nm2][j].i);
     }
  } else if (index==2) {
     idy2 = 0.5/dxy;
     for (j=oy+1;j<ny;j++) { 
        jm1=j-1;  jp1=j+1;
        for (i=ox;i<=nx;i++)  { 
           df[i][j].r = idy2*(f[i][jp1].r-f[i][jm1].r); 
           df[i][j].i = idy2*(f[i][jp1].i-f[i][jm1].i);
        }
     }
     Nm1=ny-1;  Nm2=ny-2;
     for (i=ox;i<=nx;i++)  { 
        df[i][oy].r = idy2*(-3*f[i][oy].r + 4*f[i][oy+1].r - f[i][oy+2].r);
        df[i][oy].i = idy2*(-3*f[i][oy].i + 4*f[i][oy+1].i - f[i][oy+2].i);
        df[i][ny].r = idy2*(3*f[i][ny].r - 4*f[i][Nm1].r + f[i][Nm2].r);
        df[i][ny].i = idy2*(3*f[i][ny].i - 4*f[i][Nm1].i + f[i][Nm2].i);
     }
  }
}

/* ============================================================================
   CPDIFF2D_gen:  Sept 25/93

   Find partial derivative of a complex 2D matrix with respect to index.
   using a second order finite difference scheme.  Function is not 
   necessarily defined at evenly spaced grid points.
   ============================================================================
*/
void Cpdiff2D_gen(f, x,ox,nx, y,oy,ny, index, df)
   int ox,nx,oy,ny,index;
   double *x,*y;
   dcomplex **f,**df;
{
  int i,j,Nm1,Nm2,im1,ip1,jm1,jp1;

  if (index<1 || index>2) myerror("CPDIFF2D_gen: index out of bounds");

  if (index==1) {
     Nm1=nx-1;  Nm2=nx-2;
     for (j=oy;j<=ny;j++) {
        df[ox][j].r = intplderiv(x[ox], x[ox],x[ox+1],x[ox+2],
                                      f[ox][j].r,f[ox+1][j].r,f[ox+2][j].r);   
        df[ox][j].i = intplderiv(x[ox], x[ox],x[ox+1],x[oy+2],
                                      f[ox][j].i,f[ox+1][j].i,f[ox+2][j].i);   
        df[nx][j].r = intplderiv(x[nx], x[nx],x[Nm1],x[Nm2],
                                      f[nx][j].r,f[Nm1][j].r,f[Nm2][j].r);   
        df[nx][j].i = intplderiv(x[nx], x[nx],x[Nm1],x[Nm2],
                                      f[nx][j].i,f[Nm1][j].i,f[Nm2][j].i);   
     }
     for (i=ox+1;i<nx;i++) {
        im1=i-1;  ip1=i+1;
        for (j=oy;j<=ny;j++) {
           df[i][j].r = intplderiv(x[i], x[im1],x[i],x[ip1],
                                       f[im1][j].r,f[i][j].r,f[ip1][j].r);
           df[i][j].i = intplderiv(x[i], x[im1],x[i],x[ip1],
                                       f[im1][j].i,f[i][j].i,f[ip1][j].i);
        }
     }
  } else if (index==2) {
     Nm1=ny-1;  Nm2=ny-2;
     for (i=ox;i<=nx;i++)  { 
        df[i][oy].r = intplderiv(y[oy], y[oy],y[oy+1],y[oy+2],
                                      f[i][oy].r,f[i][oy+1].r,f[i][oy+2].r);
        df[i][oy].i = intplderiv(y[oy], y[oy],y[oy+1],y[oy+2],
                                      f[i][oy].i,f[i][oy+1].i,f[i][oy+2].i);
        df[i][ny].r = intplderiv(y[ny], y[ny],y[Nm1],y[Nm2],
                                      f[i][ny].r,f[i][Nm1].r,f[i][Nm2].r);   
        df[i][ny].i = intplderiv(y[ny], y[ny],y[Nm1],y[Nm2],
                                      f[i][ny].i,f[i][Nm1].i,f[i][Nm2].i);   
     }
     for (j=oy+1;j<ny;j++) {
        jm1=j-1;  jp1=j+1;
        for (i=ox;i<=nx;i++)  { 
           df[i][j].r = intplderiv(y[j], y[jm1],y[j],y[jp1],
                                       f[i][jm1].r,f[i][j].r,f[i][jp1].r);
           df[i][j].i = intplderiv(y[j], y[jm1],y[j],y[jp1],
                                       f[i][jm1].i,f[i][j].i,f[i][jp1].i);
        }
     }
  }
}

/* =======================================================================
   CPPDIFF2D: May 1/94

   Routines to take find derivative at an interpolated point
   and to take derivative of a function fn[i](x[i]).   Method is second order 
   accurate and assumes periodicity in direction of differentiation.

   Function is assumed to be defined at evenly spaced grid points
   =======================================================================
*/
void CPpdiff2D(f,ox,nx,oy,ny, dxy,index, df)
   int ox,nx,oy,ny,index;
   double dxy;
   dcomplex **f,**df;
{
  int i,j,im1,ip1,jm1,jp1;
  double idx2,idy2;

  if (index<=0 || index>2) myerror("CPPDIFF2D: index out of bounds");

  if (index==1 && ox!=0) myerror("CPPDIFF2D: ox should equal zero");
  if (index==2 && oy!=0) myerror("CPPDIFF2D: oy should equal zero");

  if (index==1) {  /* take derivative wrt x, assume evenly spaced */
     idx2 = 0.5/dxy;
     for (i=ox+1;i<nx;i++) { 
        im1=i-1;  ip1=i+1;
        for (j=oy;j<=ny;j++) {
            df[i][j].r = idx2*(f[ip1][j].r-f[im1][j].r);
            df[i][j].i = idx2*(f[ip1][j].i-f[im1][j].i);
        }
     }
     for (j=oy;j<=ny;j++) {
        df[nx][j].r=df[ox][j].r = idx2*(f[ox+1][j].r-f[nx-1][j].r);
        df[nx][j].i=df[ox][j].i = idx2*(f[ox+1][j].i-f[nx-1][j].i);
     }
  } else if (index==2) {    /* take derivative wrt y, assume evenly spaced */
     idy2 = 0.5/dxy;
     for (j=oy+1;j<ny;j++) {
        jm1=j-1;  jp1=j+1;
        for (i=ox;i<=nx;i++) {
           df[i][j].r = idy2*(f[i][jp1].r-f[i][jm1].r);
           df[i][j].i = idy2*(f[i][jp1].i-f[i][jm1].i);
        }
     }
     for (i=ox;i<=nx;i++) {
        df[i][ny].r=df[i][oy].r = idy2*(f[i][oy+1].r-f[i][ny-1].r);
        df[i][ny].i=df[i][oy].i = idy2*(f[i][oy+1].i-f[i][ny-1].i);
     }
  }
}

/* ============================================================================
   CFPDIFF2D:  Sept 5/93

   Find partial derivative of a complex fourier 2D matrix.
   The derivative is taken of the index'th component of a spectral field.
     i.e.  diff( sum(n=0..nnx, A_n exp(I n alpha x)) ) =
              sum(n=0..nnx, I n alpha A_n exp(I n alpha x))

   ============================================================================
*/
void CFpdiff2D(f,ox,nx,oy,ny, kxy,index, df)
   int ox,nx,oy,ny,index;
   double *kxy;
   dcomplex **f,**df;
{
  int n,m;

  if (index<1 || index>2) myerror("CFPDIFF2D: index out of bounds");

  if (index==1) {    /* take spectral derivative wrt x */
     for (n=ox;n<=nx;n++)
     for (m=oy;m<=ny;m++) {
        df[n][m].r = -kxy[n]*f[n][m].i;
        df[n][m].i = kxy[n]*f[n][m].r;
     }
  } else if (index==2) {
     for (n=ox;n<=nx;n++) 
     for (m=oy;m<=ny;m++) { 
        df[n][m].r = -kxy[m]*f[n][m].i;
        df[n][m].i = kxy[m]*f[n][m].r;
     }
  }
}

/* ============================================================================
   CFPDIFF2D_gen:  Sept 5/93

   Find partial derivative of a complex fourier 2D matrix.
   The derivative is taken of the index'th component of a spectral field.
     i.e.  diff( sum(n=0..nnx, A_n exp(I n alpha x)) ) =
              sum(n=0..nnx, I n alpha A_n exp(I n alpha x))

   ============================================================================
*/
void CFpdiff2D_gen(f,ox,nx,oy,ny, kxy,index, df)
   int ox,nx,oy,ny,index;
   double *kxy;
   dcomplex **f,**df;
{
  int n,m;

  if (index<1 || index>2) myerror("CFPDIFF2D_gen: index out of bounds");

  if (index==1) {    /* take spectral derivative wrt x */
     for (n=ox;n<=nx;n++)
     for (m=oy;m<=ny;m++) {
        df[n][m].r = -kxy[n]*f[n][m].i;
        df[n][m].i = kxy[n]*f[n][m].r;
     }
  } else if (index==2) {
     for (n=ox;n<=nx;n++) 
     for (m=oy;m<=ny;m++) { 
        df[n][m].r = -kxy[m]*f[n][m].i;
        df[n][m].i = kxy[m]*f[n][m].r;
     }
  }
}

/* ============================================================================
   ZzPDIFF2DO1:  Dec 22/99

   Find partial derivative of a complex 2D matrix with respect to index.
   using a _first_ order finite difference scheme.  Derivative is presumed
   to be calculated at the intermediate gridpoint.  Direction of the
   index depends on the sign of index (i.e. see comments in PDIFF2DO1).
   By default, where df is not defined zeros are inserted.
   ============================================================================
*/
void Zzpdiff2DO1(f,ox,nx,oy,ny, dxy,index, df)
   int ox,nx,oy,ny,index;
   double dxy;
   zomplex **f,**df;
{
  int aind;
  int i,j,ip1,jp1,im1,jm1;
  double idx,idy;

  aind = iabs(index);
  if (aind==0 || aind>2) myerror("ZPDIFF2DO1: index out of bounds");

  if (index==1 && ox!=1) myerror("ZPDIFF2DO1: ox should equal one");
  if (index==2 && oy!=1) myerror("ZPDIFF2DO1: oy should equal one");
  if (index==-1 && ox!=0) myerror("ZPDIFF2DO1: ox should equal zero");
  if (index==-2 && oy!=0) myerror("ZPDIFF2DO1: oy should equal zero");

  if (index==1) {
     idx = 1.0/dxy;
     for (i=1;i<nx;i++) {
        ip1=i+1;
        for (j=oy;j<=ny;j++) {
           df[i][j].re = idx*(f[ip1][j].re-f[i][j].re);
           df[i][j].im = idx*(f[ip1][j].im-f[i][j].im);
        }
     }
     for (j=oy;j<=ny;j++) {
        df[0][j].re = df[1][j].re;
        df[0][j].im = df[1][j].im;
        df[nx][j].re = df[nx-1][j].re;
        df[nx][j].im = df[nx-1][j].im;
     }
  } else if (index==2) {
     idy = 1.0/dxy;
     for (j=1;j<ny;j++) { 
        jp1=j+1;
        for (i=ox;i<=nx;i++)  { 
           df[i][j].re = idy*(f[i][jp1].re-f[i][j].re); 
           df[i][j].im = idy*(f[i][jp1].im-f[i][j].im);
        }
     }
     for (i=ox;i<=nx;i++) {
        df[i][0].re = df[i][1].re;
        df[i][0].im = df[i][1].im;
        df[i][ny].re = df[i][ny-1].re;
        df[i][ny].im = df[i][ny-1].im;
     }
  } else if (index==-1) {
     idx = 1.0/dxy;
     for (i=1;i<=nx;i++) {
        im1=i-1;
        for (j=oy;j<=ny;j++) {
           df[i][j].re = idx*(f[i][j].re-f[im1][j].re);
           df[i][j].im = idx*(f[i][j].im-f[im1][j].im);
        }
     }
  } else if (index==-2) {
     idy = 1.0/dxy;
     for (j=1;j<=ny;j++) { 
        jm1=j-1;
        for (i=ox;i<=nx;i++)  { 
           df[i][j].re = idy*(f[i][j].re-f[i][jm1].re); 
           df[i][j].im = idy*(f[i][j].im-f[i][jm1].im);
        }
     }
  }
}

/* ============================================================================
   ZPPDIFF2DO1:  Dec 7/93

   Find partial derivative of a complex 2D matrix with respect to index.
   using a _first_ order finite difference scheme.  Derivative is presumed
   to be calculated at the intermediate gridpoint.  Direction of the
   index depends on the sign of index (i.e. see comments in PDIFF2DO1).
   By default, where df is not defined zeros are inserted.

   This procedure differs from ZPDIFF2DO1 in that the domain is assumed to
   be periodic at the boundary.
   ============================================================================
*/
void ZPpdiff2DO1(f,ox,nx,oy,ny, dxy,index, df)
   int ox,nx,oy,ny,index;
   double dxy;
   zomplex **f,**df;
{
  int aind;
  int i,j,ip1,jp1,im1,jm1;
  double idx,idy;

  aind = iabs(index);
  if (aind==0 || aind>2) myerror("ZPPDIFF2DO1: index out of bounds");

  if (aind==1 && ox!=0) myerror("ZPPDIFF2DO1: ox should equal zero");
  if (aind==2 && oy!=0) myerror("ZPPDIFF2DO1: oy should equal zero");

  if (index==1) {
     idx = 1.0/dxy;
     for (i=0;i<nx;i++) {
        ip1=i+1;
        for (j=oy;j<=ny;j++) {
           df[i][j].re = idx*(f[ip1][j].re-f[i][j].re);
           df[i][j].im = idx*(f[ip1][j].im-f[i][j].im);
        }
     }
     for (j=oy;j<=ny;j++) df[nx][j] = df[0][j];
  } else if (index==2) {
     idy = 1.0/dxy;
     for (j=0;j<ny;j++) { 
        jp1=j+1;
        for (i=ox;i<=nx;i++)  { 
           df[i][j].re = idy*(f[i][jp1].re-f[i][j].re); 
           df[i][j].im = idy*(f[i][jp1].im-f[i][j].im);
        }
     }
     for (i=ox;i<=nx;i++) df[i][ny] = df[i][0];
  } else if (index==-1) {
     idx = 1.0/dxy;
     for (i=1;i<=nx;i++) {
        im1=i-1;
        for (j=oy;j<=ny;j++) {
           df[i][j].re = idx*(f[i][j].re-f[im1][j].re);
           df[i][j].im = idx*(f[i][j].im-f[im1][j].im);
        }
     }
     for (j=oy;j<=ny;j++) df[0][j] = df[nx][j];
  } else if (index==-2) {
     idy = 1.0/dxy;
     for (j=1;j<=ny;j++) { 
        jm1=j-1;
        for (i=ox;i<=nx;i++)  { 
           df[i][j].re = idy*(f[i][j].re-f[i][jm1].re); 
           df[i][j].im = idy*(f[i][j].im-f[i][jm1].im);
        }
     }
     for (i=ox;i<=nx;i++) df[i][0] = df[i][ny];
  }
}

/* ============================================================================
   ZPDIFF2D:  Sept 5/93

   Find partial derivative of a complex 2D matrix with respect to index.
   using a second order finite difference scheme.  Function is defined
   at evenly spaced grid points with spacing dxy in given direction of
   differentiation.
   ============================================================================
*/
void Zpdiff2D(f,ox,nx,oy,ny, dxy,index, df)
   int ox,nx,oy,ny,index;
   double dxy;
   zomplex **f,**df;
{
  int i,j,Nm1,Nm2,im1,ip1,jm1,jp1;
  double idx2,idy2;

  if (index<1 || index>2) myerror("ZPDIFF2D: index out of bounds");

  if (index==1) {
     idx2 = 0.5/dxy;
     for (i=ox+1;i<nx;i++) {
        im1=i-1;  ip1=i+1;
        for (j=oy;j<=ny;j++) {
           df[i][j].re = idx2*(f[ip1][j].re-f[im1][j].re);
           df[i][j].im = idx2*(f[ip1][j].im-f[im1][j].im);
        }
     }
     Nm1=nx-1;  Nm2=nx-2;
     for (j=oy;j<=ny;j++) {
        df[ox][j].re = idx2*(-3*f[ox][j].re + 4*f[ox+1][j].re - f[ox+2][j].re);
        df[ox][j].im = idx2*(-3*f[ox][j].im + 4*f[ox+1][j].im - f[ox+2][j].im);
        df[nx][j].re = idx2*(3*f[nx][j].re - 4*f[Nm1][j].re + f[Nm2][j].re);
        df[nx][j].im = idx2*(3*f[nx][j].im - 4*f[Nm1][j].im + f[Nm2][j].im);
     }
  } else if (index==2) {
     idy2 = 0.5/dxy;
     for (j=oy+1;j<ny;j++) { 
        jm1=j-1;  jp1=j+1;
        for (i=ox;i<=nx;i++)  { 
           df[i][j].re = idy2*(f[i][jp1].re-f[i][jm1].re); 
           df[i][j].im = idy2*(f[i][jp1].im-f[i][jm1].im);
        }
     }
     Nm1=ny-1;  Nm2=ny-2;
     for (i=ox;i<=nx;i++)  { 
        df[i][oy].re = idy2*(-3*f[i][oy].re + 4*f[i][oy+1].re - f[i][oy+2].re);
        df[i][oy].im = idy2*(-3*f[i][oy].im + 4*f[i][oy+1].im - f[i][oy+2].im);
        df[i][ny].re = idy2*(3*f[i][ny].re - 4*f[i][Nm1].re + f[i][Nm2].re);
        df[i][ny].im = idy2*(3*f[i][ny].im - 4*f[i][Nm1].im + f[i][Nm2].im);
     }
  }
}

/* ============================================================================
   ZPDIFF2D_gen:  Sept 25/93

   Find partial derivative of a complex 2D matrix with respect to index.
   using a second order finite difference scheme.  Function is not 
   necessarily defined at evenly spaced grid points.
   ============================================================================
*/
void Zpdiff2D_gen(f, x,ox,nx, y,oy,ny, index, df)
   int ox,nx,oy,ny,index;
   double *x,*y;
   zomplex **f,**df;
{
  int i,j,Nm1,Nm2,im1,ip1,jm1,jp1;

  if (index<1 || index>2) myerror("ZPDIFF2D_gen: index out of bounds");

  if (index==1) { /* x-derivative */
     Nm1=nx-1;  Nm2=nx-2;
     for (j=oy;j<=ny;j++) {
        df[ox][j].re = intplderiv(x[ox], x[ox],x[ox+1],x[ox+2],
                                      f[ox][j].re,f[ox+1][j].re,f[ox+2][j].re);   
        df[ox][j].im = intplderiv(x[ox], x[ox],x[ox+1],x[oy+2],
                                      f[ox][j].im,f[ox+1][j].im,f[ox+2][j].im);   
        df[nx][j].re = intplderiv(x[nx], x[nx],x[Nm1],x[Nm2],
                                      f[nx][j].re,f[Nm1][j].re,f[Nm2][j].re);   
        df[nx][j].im = intplderiv(x[nx], x[nx],x[Nm1],x[Nm2],
                                      f[nx][j].im,f[Nm1][j].im,f[Nm2][j].im);   
     }
     for (i=ox+1;i<nx;i++) {
        im1=i-1;  ip1=i+1;
        for (j=oy;j<=ny;j++) {
           df[i][j].re = intplderiv(x[i], x[im1],x[i],x[ip1],
                                       f[im1][j].re,f[i][j].re,f[ip1][j].re);
           df[i][j].im = intplderiv(x[i], x[im1],x[i],x[ip1],
                                       f[im1][j].im,f[i][j].im,f[ip1][j].im);
        }
     }
  } else if (index==2) { /* y/z-derivative */
     Nm1=ny-1;  Nm2=ny-2;
     for (i=ox;i<=nx;i++)  { 
        df[i][oy].re = intplderiv(y[oy], y[oy],y[oy+1],y[oy+2],
                                      f[i][oy].re,f[i][oy+1].re,f[i][oy+2].re);
        df[i][oy].im = intplderiv(y[oy], y[oy],y[oy+1],y[oy+2],
                                      f[i][oy].im,f[i][oy+1].im,f[i][oy+2].im);
        df[i][ny].re = intplderiv(y[ny], y[ny],y[Nm1],y[Nm2],
                                      f[i][ny].re,f[i][Nm1].re,f[i][Nm2].re);   
        df[i][ny].im = intplderiv(y[ny], y[ny],y[Nm1],y[Nm2],
                                      f[i][ny].im,f[i][Nm1].im,f[i][Nm2].im);   
     }
     for (j=oy+1;j<ny;j++) {
        jm1=j-1;  jp1=j+1;
        for (i=ox;i<=nx;i++)  { 
           df[i][j].re = intplderiv(y[j], y[jm1],y[j],y[jp1],
                                       f[i][jm1].re,f[i][j].re,f[i][jp1].re);
           df[i][j].im = intplderiv(y[j], y[jm1],y[j],y[jp1],
                                       f[i][jm1].im,f[i][j].im,f[i][jp1].im);
        }
     }
  }
}

/* =======================================================================
   ZPPDIFF2D: May 1/94

   Routines to take find derivative at an interpolated point
   and to take derivative of a function fn[i](x[i]).   Method is second order 
   accurate and assumes periodicity in direction of differentiation.

   Function is assumed to be defined at evenly spaced grid points
   =======================================================================
*/
void ZPpdiff2D(f,ox,nx,oy,ny, dxy,index, df)
   int ox,nx,oy,ny,index;
   double dxy;
   zomplex **f,**df;
{
  int i,j,im1,ip1,jm1,jp1;
  double idx2,idy2;

  if (index<=0 || index>2) myerror("ZPPDIFF2D: index out of bounds");

  if (index==1 && ox!=0) myerror("ZPPDIFF2D: ox should equal zero");
  if (index==2 && oy!=0) myerror("ZPPDIFF2D: oy should equal zero");

  if (index==1) {  /* take derivative wrt x, assume evenly spaced */
     idx2 = 0.5/dxy;
     for (i=ox+1;i<nx;i++) { 
        im1=i-1;  ip1=i+1;
        for (j=oy;j<=ny;j++) {
            df[i][j].re = idx2*(f[ip1][j].re-f[im1][j].re);
            df[i][j].im = idx2*(f[ip1][j].im-f[im1][j].im);
        }
     }
     for (j=oy;j<=ny;j++) {
        df[nx][j].re=df[ox][j].re = idx2*(f[ox+1][j].re-f[nx-1][j].re);
        df[nx][j].im=df[ox][j].im = idx2*(f[ox+1][j].im-f[nx-1][j].im);
     }
  } else if (index==2) {    /* take derivative wrt y, assume evenly spaced */
     idy2 = 0.5/dxy;
     for (j=oy+1;j<ny;j++) {
        jm1=j-1;  jp1=j+1;
        for (i=ox;i<=nx;i++) {
           df[i][j].re = idy2*(f[i][jp1].re-f[i][jm1].re);
           df[i][j].im = idy2*(f[i][jp1].im-f[i][jm1].im);
        }
     }
     for (i=ox;i<=nx;i++) {
        df[i][ny].re=df[i][oy].re = idy2*(f[i][oy+1].re-f[i][ny-1].re);
        df[i][ny].im=df[i][oy].im = idy2*(f[i][oy+1].im-f[i][ny-1].im);
     }
  }
}

/* ============================================================================
   ZFPDIFF2D:  Sept 5/93

   Find partial derivative of a complex fourier 2D matrix.
   The derivative is taken of the index'th component of a spectral field.
     i.e.  diff( sum(n=0..nnx, A_n exp(I n alpha x)) ) =
              sum(n=0..nnx, I n alpha A_n exp(I n alpha x))

   ============================================================================
*/
void ZFpdiff2D(f,ox,nx,oy,ny, kxy,index, df)
   int ox,nx,oy,ny,index;
   double *kxy;
   zomplex **f,**df;
{
  int n,m;

  if (index<1 || index>2) myerror("ZFPDIFF2D: index out of bounds");

  if (index==1) {    /* take spectral derivative wrt x */
     for (n=ox;n<=nx;n++)
     for (m=oy;m<=ny;m++) {
        df[n][m].re = -kxy[n]*f[n][m].im;
        df[n][m].im = kxy[n]*f[n][m].re;
     }
  } else if (index==2) {
     for (n=ox;n<=nx;n++) 
     for (m=oy;m<=ny;m++) { 
        df[n][m].re = -kxy[m]*f[n][m].im;
        df[n][m].im = kxy[m]*f[n][m].re;
     }
  }
}

/* ============================================================================
   ZFPDIFF2D_gen:  Sept 5/93

   Find partial derivative of a complex fourier 2D matrix.
   The derivative is taken of the index'th component of a spectral field.
     i.e.  diff( sum(n=0..nnx, A_n exp(I n alpha x)) ) =
              sum(n=0..nnx, I n alpha A_n exp(I n alpha x))

   ============================================================================
*/
void ZFpdiff2D_gen(f,ox,nx,oy,ny, kxy,index, df)
   int ox,nx,oy,ny,index;
   double *kxy;
   zomplex **f,**df;
{
  int n,m;

  if (index<1 || index>2) myerror("ZFPDIFF2D_gen: index out of bounds");

  if (index==1) {    /* take spectral derivative wrt x */
     for (n=ox;n<=nx;n++)
     for (m=oy;m<=ny;m++) {
        df[n][m].re = -kxy[n]*f[n][m].im;
        df[n][m].im = kxy[n]*f[n][m].re;
     }
  } else if (index==2) {
     for (n=ox;n<=nx;n++) 
     for (m=oy;m<=ny;m++) { 
        df[n][m].re = -kxy[m]*f[n][m].im;
        df[n][m].im = kxy[m]*f[n][m].re;
     }
  }
}

/* ============================================================================
   PDIFF3DO1:  Dec 15/93

   Find partial derivative with respect to index of a real field defined 
   at evenly spaced grid points on a 3D matrix using a _first_ order finite 
   difference scheme.  Derivative is presumed to be calculated at the 
   intermediate gridpoint.   The direction of differentiation depends on
   the sign of index (i.e. see comments to PDIFF2DO1).
   By default, zeros are inserted in df where it is not defined.
   ============================================================================
*/
void pdiff3DO1(f,ox,nx,oy,ny,oz,nz, dxyz,index, df)
   int ox,nx,oy,ny,oz,nz,index;
   double dxyz;
   double ***f,***df;
{
  int aind;
  int i,j,k;
  int ip1,jp1,kp1,im1,jm1,km1;
  double idx,idy,idz;

  aind = iabs(index);
  if (aind==0 || aind>3) myerror("RPDIFF3DO1: index out of bounds");

  if (index==1 && ox!=1) myerror("RPDIFF3DO1: ox should equal one");
  if (index==2 && oy!=1) myerror("RPDIFF3DO1: oy should equal one");
  if (index==3 && oz!=1) myerror("RPDIFF3DO1: oz should equal one");
  if (index==-1 && ox!=0) myerror("RPDIFF3DO1: ox should equal zero");
  if (index==-2 && oy!=0) myerror("RPDIFF3DO1: oy should equal zero");
  if (index==-3 && oz!=0) myerror("RPDIFF3DO1: oz should equal zero");

  if (index==1) {
     idx = 1.0/dxyz;
     for (i=1;i<nx;i++) {
        ip1=i+1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) df[i][j][k] = idx*(f[ip1][j][k]-f[i][j][k]); 
     }
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) {
        df[0][j][k] = df[1][j][k];
        df[nx][j][k] = df[nx-1][j][k];
     }
  } else if (index==2) {
     idy = 1.0/dxyz;
     for (j=1;j<ny;j++) {
        jp1=j+1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k] = idy*(f[i][jp1][k]-f[i][j][k]); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) {
        df[i][0][k] = df[i][1][k];
        df[i][ny][k] = df[i][ny-1][k];
     }
  } else if (index==3) {
     idz = 1.0/dxyz;
     for (k=1;k<nz;k++) {
        kp1=k+1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k] = idz*(f[i][j][kp1]-f[i][j][k]); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (j=oy;j<=ny;j++) {
        df[i][j][0] = df[i][j][1];
        df[i][j][nz] = df[i][j][nz-1];
     }
  } else if (index==-1) {
     idx = 1.0/dxyz;
     for (i=1;i<=nx;i++) {
        im1=i-1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k] = idx*(f[i][j][k]-f[im1][j][k]); 
        }
     }
  } else if (index==-2) {
     idy = 1.0/dxyz;
     for (j=1;j<=ny;j++) {
        jm1=j-1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k] = idy*(f[i][j][k]-f[i][jm1][k]); 
        }
     }
  } else if (index==-3) {
     idz = 1.0/dxyz;
     for (k=1;k<=nz;k++) {
        km1=k-1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k] = idz*(f[i][j][k]-f[i][j][km1]); 
        }
     }
  }
}

/* ============================================================================
   PPDIFF3DO1:  Dec 15/93

   Find partial derivative with respect to index of a real field defined 
   at evenly spaced grid points on a 3D matrix using a _first_ order finite 
   difference scheme.  Derivative is presumed to be calculated at the 
   intermediate gridpoint.   The direction of differentiation depends on
   the sign of index (i.e. see comments to PDIFF2DO1).
   By default, zeros are inserted in df where it is not defined.

   This procedure differs from PDIFFO1 in that is assumes the domain
   has periodic boundaries.
   ============================================================================
*/
void Ppdiff3DO1(f,ox,nx,oy,ny,oz,nz, dxyz,index, df)
   int ox,nx,oy,ny,oz,nz,index;
   double dxyz;
   double ***f,***df;
{
  int aind;
  int i,j,k;
  int ip1,jp1,kp1,im1,jm1,km1;
  double idx,idy,idz;

  aind = iabs(index);
  if (aind==0 || aind>3) myerror("RPPDIFF3DO1: index out of bounds");

  if (aind==1 && ox!=0) myerror("RPPDIFF3DO1: ox should equal zero");
  if (aind==2 && oy!=0) myerror("RPPDIFF3DO1: oy should equal zero");
  if (aind==3 && oz!=0) myerror("RPPDIFF3DO1: oz should equal zero");

  if (index==1) {
     idx = 1.0/dxyz;
     for (i=0;i<nx;i++) {
        ip1=i+1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k] = idx*(f[ip1][j][k]-f[i][j][k]); 
        }
     }
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) df[nx][j][k] = df[0][j][k];
  } else if (index==2) {
     idy = 1.0/dxyz;
     for (j=0;j<ny;j++) {
        jp1=j+1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k] = idy*(f[i][jp1][k]-f[i][j][k]); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) df[i][ny][k] = df[i][0][k];
  } else if (index==3) {
     idz = 1.0/dxyz;
     for (k=0;k<nz;k++) {
        kp1=k+1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k] = idz*(f[i][j][kp1]-f[i][j][k]); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (j=oy;j<=ny;j++) df[i][j][nz] = df[i][j][0];
  } else if (index==-1) {
     idx = 1.0/dxyz;
     for (i=1;i<=nx;i++) {
        im1=i-1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k] = idx*(f[i][j][k]-f[im1][j][k]); 
        }
     }
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) df[0][j][k] = df[nx][j][k];
  } else if (index==-2) {
     idy = 1.0/dxyz;
     for (j=1;j<=ny;j++) {
        jm1=j-1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k] = idy*(f[i][j][k]-f[i][jm1][k]); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) df[i][0][k] = df[i][ny][k];
  } else if (index==-3) {
     idz = 1.0/dxyz;
     for (k=1;k<=nz;k++) {
        km1=k-1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k] = idz*(f[i][j][k]-f[i][j][km1]); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (j=oy;j<=ny;j++) df[i][j][0] = df[i][j][nz];
  }
}

/* ============================================================================
   ZPDIFF3DO1:  Feb 16/94

   Find partial derivative with respect to index of a real field defined 
   at evenly spaced grid points on a 3D matrix using a _first_ order finite 
   difference scheme.  Derivative is presumed to be calculated at the 
   intermediate gridpoint.   The direction of differentiation depends on
   the sign of index (i.e. see comments to PDIFF2DO1).
   By default, zeros are inserted in df where it is not defined.

   This procedure differs from PDIFFO1 in that is assumes the domain
   has periodic boundaries.
   ============================================================================
*/
void Zpdiff3DO1(f,ox,nx,oy,ny,oz,nz, dxyz,index, df)
   int ox,nx,oy,ny,oz,nz,index;
   double dxyz;
   double ***f,***df;
{
  int aind;
  int i,j,k;
  int ip1,jp1,kp1,im1,jm1,km1;
  double idx,idy,idz;

  aind = iabs(index);
  if (aind==0 || aind>3) myerror("ZPDIFF3DO1: index out of bounds");

  if (index==1 && ox!=1) myerror("ZPDIFF3DO1: ox should equal one");
  if (index==2 && oy!=1) myerror("ZPDIFF3DO1: oy should equal one");
  if (index==3 && oz!=1) myerror("ZPDIFF3DO1: oz should equal one");
  if (index==-1 && ox!=0) myerror("ZPDIFF3DO1: ox should equal zero");
  if (index==-2 && oy!=0) myerror("ZPDIFF3DO1: oy should equal zero");
  if (index==-3 && oz!=0) myerror("ZPDIFF3DO1: oz should equal zero");

  if (index==1) {
     idx = 1.0/dxyz;
     for (i=0;i<nx;i++) {
        ip1=i+1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k] = idx*(f[ip1][j][k]-f[i][j][k]); 
        }
     }
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) df[nx][j][k] = df[0][j][k] = 0.0;
  } else if (index==2) {
     idy = 1.0/dxyz;
     for (j=0;j<ny;j++) {
        jp1=j+1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k] = idy*(f[i][jp1][k]-f[i][j][k]); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) df[i][ny][k] = df[i][0][k] = 0.0;
  } else if (index==3) {
     idz = 1.0/dxyz;
     for (k=0;k<nz;k++) {
        kp1=k+1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k] = idz*(f[i][j][kp1]-f[i][j][k]); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (j=oy;j<=ny;j++) df[i][j][nz] = df[i][j][0] = 0.0;
  } else if (index==-1) {
     idx = 1.0/dxyz;
     for (i=1;i<=nx;i++) {
        im1=i-1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k] = idx*(f[i][j][k]-f[im1][j][k]); 
        }
     }
  } else if (index==-2) {
     idy = 1.0/dxyz;
     for (j=1;j<=ny;j++) {
        jm1=j-1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k] = idy*(f[i][j][k]-f[i][jm1][k]); 
        }
     }
  } else if (index==-3) {
     idz = 1.0/dxyz;
     for (k=1;k<=nz;k++) {
        km1=k-1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k] = idz*(f[i][j][k]-f[i][j][km1]); 
        }
     }
  }
}

/* ============================================================================
   CPDIFF3DO1:  Sept 20/93

   Find partial derivative with respect to index of a complex field defined 
   at evenly spaced grid points on a 3D matrix using a _first_ order finite 
   difference scheme.  Derivative is presumed to be calculated at the 
   intermediate gridpoint.   The direction of differentiation depends on
   the sign of index (i.e. see comments to PDIFF2DO1).
   By default, zeros are inserted in df where it is not defined.
   ============================================================================
*/
void Cpdiff3DO1(f,ox,nx,oy,ny,oz,nz, dxyz,index, df)
   int ox,nx,oy,ny,oz,nz,index;
   double dxyz;
   dcomplex ***f,***df;
{
  int aind;
  int i,j,k;
  int ip1,jp1,kp1,im1,jm1,km1;
  double idx,idy,idz;

  aind = iabs(index);
  if (aind==0 || aind>3) myerror("CPDIFF3DO1: index out of bounds");

  if (index==1 && ox!=1) myerror("CPDIFF3DO1: ox should equal one");
  if (index==2 && oy!=1) myerror("CPDIFF3DO1: oy should equal one");
  if (index==3 && oz!=1) myerror("CPDIFF3DO1: oz should equal one");
  if (index==-1 && ox!=0) myerror("CPDIFF3DO1: ox should equal zero");
  if (index==-2 && oy!=0) myerror("CPDIFF3DO1: oy should equal zero");
  if (index==-3 && oz!=0) myerror("CPDIFF3DO1: oz should equal zero");

  if (index==1) {
     idx = 1.0/dxyz;
     for (i=1;i<nx;i++) {
        ip1=i+1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].r = idx*(f[ip1][j][k].r-f[i][j][k].r); 
           df[i][j][k].i = idx*(f[ip1][j][k].i-f[i][j][k].i); 
        }
     }
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) {
        df[0][j][k].r = df[1][j][k].r;
        df[0][j][k].i = df[1][j][k].i;
        df[nx][j][k].r = df[nx-1][j][k].r;
        df[nx][j][k].i = df[nx-1][j][k].i;
     }
  } else if (index==2) {
     idy = 1.0/dxyz;
     for (j=1;j<ny;j++) {
        jp1=j+1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].r = idy*(f[i][jp1][k].r-f[i][j][k].r); 
           df[i][j][k].i = idy*(f[i][jp1][k].i-f[i][j][k].i); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) {
        df[i][0][k].r = df[i][1][k].r;
        df[i][0][k].i = df[i][1][k].i;
        df[i][ny][k].r = df[i][ny-1][k].r;
        df[i][ny][k].i = df[i][ny-1][k].i;
     }
  } else if (index==3) {
     idz = 1.0/dxyz;
     for (k=1;k<nz;k++) {
        kp1=k+1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k].r = idz*(f[i][j][kp1].r-f[i][j][k].r); 
           df[i][j][k].i = idz*(f[i][j][kp1].i-f[i][j][k].i); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (j=oy;j<=ny;j++) {
        df[i][j][0].r = df[i][j][1].r;
        df[i][j][0].i = df[i][j][1].i;
        df[i][j][nz].r = df[i][j][nz-1].r;
        df[i][j][nz].i = df[i][j][nz-1].i;
     }
  } else if (index==-1) {
     idx = 1.0/dxyz;
     for (i=1;i<=nx;i++) {
        im1=i-1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].r = idx*(f[i][j][k].r-f[im1][j][k].r); 
           df[i][j][k].i = idx*(f[i][j][k].i-f[im1][j][k].i); 
        }
     }
  } else if (index==-2) {
     idy = 1.0/dxyz;
     for (j=1;j<=ny;j++) {
        jm1=j-1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].r = idy*(f[i][j][k].r-f[i][jm1][k].r); 
           df[i][j][k].i = idy*(f[i][j][k].i-f[i][jm1][k].i); 
        }
     }
  } else if (index==-3) {
     idz = 1.0/dxyz;
     for (k=1;k<=nz;k++) {
        km1=k-1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k].r = idz*(f[i][j][k].r-f[i][j][km1].r); 
           df[i][j][k].i = idz*(f[i][j][k].i-f[i][j][km1].i); 
        }
     }
  }
}

/* ============================================================================
   CPPDIFF3DO1:  Dec 7/93

   Find partial derivative with respect to index of a complex field defined 
   at evenly spaced grid points on a 3D matrix using a _first_ order finite 
   difference scheme.  Derivative is presumed to be calculated at the 
   intermediate gridpoint.   The direction of differentiation depends on
   the sign of index (i.e. see comments to PDIFF2DO1).
   By default, zeros are inserted in df where it is not defined.

   This procedure differs from CPDIFFO1 in that is assumes the domain
   has periodic boundaries.
   ============================================================================
*/

void CPpdiff3DO1(f,ox,nx,oy,ny,oz,nz, dxyz,index, df)
   int ox,nx,oy,ny,oz,nz,index;
   double dxyz;
   dcomplex ***f,***df;
{
  int aind;
  int i,j,k;
  int ip1,jp1,kp1,im1,jm1,km1;
  double idx,idy,idz;

  aind = iabs(index);
  if (aind==0 || aind>3) myerror("CPPDIFF3DO1: index out of bounds");

  if (aind==1 && ox!=0) myerror("CPPDIFF3DO1: ox should equal zero");
  if (aind==2 && oy!=0) myerror("CPPDIFF3DO1: oy should equal zero");
  if (aind==3 && oz!=0) myerror("CPPDIFF3DO1: oz should equal zero");

  if (index==1) {
     idx = 1.0/dxyz;
     for (i=0;i<nx;i++) {
        ip1=i+1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].r = idx*(f[ip1][j][k].r-f[i][j][k].r); 
           df[i][j][k].i = idx*(f[ip1][j][k].i-f[i][j][k].i); 
        }
     }
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) df[nx][j][k] = df[0][j][k];
  } else if (index==2) {
     idy = 1.0/dxyz;
     for (j=0;j<ny;j++) {
        jp1=j+1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].r = idy*(f[i][jp1][k].r-f[i][j][k].r); 
           df[i][j][k].i = idy*(f[i][jp1][k].i-f[i][j][k].i); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) df[i][ny][k] = df[i][0][k];
  } else if (index==3) {
     idz = 1.0/dxyz;
     for (k=0;k<nz;k++) {
        kp1=k+1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k].r = idz*(f[i][j][kp1].r-f[i][j][k].r); 
           df[i][j][k].i = idz*(f[i][j][kp1].i-f[i][j][k].i); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (j=oy;j<=ny;j++) df[i][j][nz] = df[i][j][0];
  } else if (index==-1) {
     idx = 1.0/dxyz;
     for (i=1;i<=nx;i++) {
        im1=i-1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].r = idx*(f[i][j][k].r-f[im1][j][k].r); 
           df[i][j][k].i = idx*(f[i][j][k].i-f[im1][j][k].i); 
        }
     }
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) df[0][j][k] = df[nx][j][k];
  } else if (index==-2) {
     idy = 1.0/dxyz;
     for (j=1;j<=ny;j++) {
        jm1=j-1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].r = idy*(f[i][j][k].r-f[i][jm1][k].r); 
           df[i][j][k].i = idy*(f[i][j][k].i-f[i][jm1][k].i); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) df[i][0][k] = df[i][ny][k];
  } else if (index==-3) {
     idz = 1.0/dxyz;
     for (k=1;k<=nz;k++) {
        km1=k-1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k].r = idz*(f[i][j][k].r-f[i][j][km1].r); 
           df[i][j][k].i = idz*(f[i][j][k].i-f[i][j][km1].i); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (j=oy;j<=ny;j++) df[i][j][0] = df[i][j][nz];
  }
}

/* ============================================================================
   CPDIFF3D:  Sept 20/93

   Find partial derivative of a complex 3D matrix with respect to index
   using a second order finite difference scheme.  Function is defined
   on an evenly spaced grid.
   ============================================================================
*/

void Cpdiff3D(f,ox,nx,oy,ny,oz,nz, dxyz,index, df)
   int ox,nx,oy,ny,oz,nz,index;
   double dxyz;
   dcomplex ***f,***df;
{
  int i,j,k;
  int ip1,im1,jp1,jm1,kp1,km1,Nm1,Nm2;
  double idx2,idy2,idz2;

  if (index<1 || index>3) myerror("CPDIFF3D: index out of bounds");

  if (index==1) {
     idx2 = 0.5/dxyz;
     Nm1=nx-1;   Nm2=nx-2;
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) { 
        df[ox][j][k].r = 
           idx2*(-3*f[ox][j][k].r + 4*f[ox+1][j][k].r - f[ox+2][j][k].r);
        df[ox][j][k].i = 
           idx2*(-3*f[ox][j][k].i + 4*f[ox+1][j][k].i - f[ox+2][j][k].i);
        df[nx][j][k].r = 
           idx2*(3*f[nx][j][k].r - 4*f[Nm1][j][k].r + f[Nm2][j][k].r);
        df[nx][j][k].i = 
           idx2*(3*f[nx][j][k].i - 4*f[Nm1][j][k].i + f[Nm2][j][k].i);
     }
     for (i=ox+1;i<nx;i++) {
        ip1=i+1;  im1=i-1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].r = idx2*(f[ip1][j][k].r-f[im1][j][k].r); 
           df[i][j][k].i = idx2*(f[ip1][j][k].i-f[im1][j][k].i); 
        }
     }
  } else if (index==2) {
     idy2 = 0.5/dxyz;
     Nm1=ny-1;   Nm2=ny-2;
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) { 
        df[i][oy][k].r = 
           idy2*(-3*f[i][oy][k].r + 4*f[i][oy+1][k].r - f[i][oy+2][k].r);
        df[i][oy][k].i = 
           idy2*(-3*f[i][oy][k].i + 4*f[i][oy+1][k].i - f[i][oy+2][k].i);
        df[i][ny][k].r = 
           idy2*(3*f[i][ny][k].r - 4*f[i][Nm1][k].r + f[i][Nm2][k].r);
        df[i][ny][k].i = 
           idy2*(3*f[i][ny][k].i - 4*f[i][Nm1][k].i + f[i][Nm2][k].i);
     }
     for (j=oy+1;j<ny;j++) {
        jp1=j+1;  jm1=j-1;
        for (i=ox;i<=nx;i++)  
        for (k=ox;k<=nz;k++) { 
           df[i][j][k].r = idy2*(f[i][jp1][k].r-f[i][jm1][k].r); 
           df[i][j][k].i = idy2*(f[i][jp1][k].i-f[i][jm1][k].i); 
        }
     }
  } else if (index==3) {
     idz2 = 0.5/dxyz;
     Nm1=nz-1;   Nm2=nz-2;
     for (i=ox;i<=nx;i++)  
     for (j=oy;j<=ny;j++) {  
        df[i][j][oz].r = 
           idz2*(-3*f[i][j][oz].r + 4*f[i][j][oz+1].r - f[i][j][oz+2].r);
        df[i][j][oz].i = 
           idz2*(-3*f[i][j][oz].i + 4*f[i][j][oz+1].i - f[i][j][oz+2].i);
        df[i][j][nz].r = 
           idz2*(3*f[i][j][nz].r - 4*f[i][j][Nm1].r + f[i][j][Nm2].r);
        df[i][j][nz].i = 
           idz2*(3*f[i][j][nz].i - 4*f[i][j][Nm1].i + f[i][j][Nm2].i);
     }
     for (k=oz+1;k<nz;k++) {
        kp1=k+1; km1=k-1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k].r = idz2*(f[i][j][kp1].r-f[i][j][km1].r); 
           df[i][j][k].i = idz2*(f[i][j][kp1].i-f[i][j][km1].i); 
        }
     }
  }
}

/* ============================================================================
   CPDIFF3D_gen:  Sept 20/93

   Find partial derivative of a complex 3D matrix with respect to index
   using a second order finite difference scheme.  Function is not
   necessarily defined at evenly space points.
   ============================================================================
*/

void Cpdiff3D_gen(f, x,ox,nx, y,oy,ny, z,oz,nz, index, df)
   int ox,nx,oy,ny,oz,nz,index;
   double *x,*y,*z;
   dcomplex ***f,***df;
{
  int i,j,k;
  int ip1,im1,jp1,jm1,kp1,km1,Nm1,Nm2;

  if (index<1 || index>3) myerror("CPDIFF3D_gen: index out of bounds");

  if (index==1) {
     Nm1=nx-1; Nm2=nx-2;
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) { 
        df[ox][j][k].r = intplderiv(x[ox], 
           x[ox],x[ox+1],x[ox+2],f[ox][j][k].r,f[ox+1][j][k].r,f[ox+2][j][k].r);
        df[ox][j][k].i = intplderiv(x[ox], 
           x[ox],x[ox+1],x[ox+2],f[ox][j][k].i,f[ox+1][j][k].i,f[ox+2][j][k].i);
        df[nx][j][k].r = intplderiv(x[nx], 
           x[nx],x[Nm1],x[Nm2],f[nx][j][k].r,f[Nm1][j][k].r,f[Nm2][j][k].r);
        df[nx][j][k].i = intplderiv(x[nx],
           x[nx],x[Nm1],x[Nm2],f[nx][j][k].i,f[Nm1][j][k].i,f[Nm2][j][k].i);
     }

     for (i=ox+1;i<nx;i++) {  
        im1=i-1; ip1=i+1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].r = intplderiv(x[i], 
              x[im1],x[i],x[ip1],f[im1][j][k].r,f[i][j][k].r,f[ip1][j][k].r);
           df[i][j][k].i = intplderiv(x[i], 
              x[im1],x[i],x[ip1],f[im1][j][k].i,f[i][j][k].i,f[ip1][j][k].i);
        }
     }

  } else if (index==2) {
     Nm1=ny-1; Nm2=ny-2;
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) { 
        df[i][oy][k].r = intplderiv(y[oy], 
           y[oy],y[oy+1],y[oy+2],f[i][oy][k].r,f[i][oy+1][k].r,f[i][oy+2][k].r);
        df[i][oy][k].i = intplderiv(y[oy], 
           y[oy],y[oy+1],y[oy+2],f[i][oy][k].i,f[i][oy+1][k].i,f[i][oy+2][k].i);
        df[i][ny][k].r = intplderiv(y[ny], 
           y[ny],y[Nm1],y[Nm2],f[i][ny][k].r,f[i][Nm1][k].r,f[i][Nm2][k].r);
        df[i][ny][k].i = intplderiv(y[ny],
           y[ny],y[Nm1],y[Nm2],f[i][ny][k].i,f[i][Nm1][k].i,f[i][Nm2][k].i);
     }

     for (j=oy+1;j<ny;j++) {  
        jm1=j-1;  jp1=j+1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].r = intplderiv(y[j], 
              y[jm1],y[j],y[jp1],f[i][jm1][k].r,f[i][j][k].r,f[i][jp1][k].r);
           df[i][j][k].i = intplderiv(y[j], 
              y[jm1],y[j],y[jp1],f[i][jm1][k].i,f[i][j][k].i,f[i][jp1][k].i);
        }
     }

  } else if (index==3) {
     Nm1=nz-1;  Nm2=nz-2;
     for (i=ox;i<=nx;i++)
     for (j=oy;j<=ny;j++) {   
        df[i][j][oz].r = intplderiv(z[oz], 
           z[oz],z[oz+1],z[oz+2],f[i][j][oz].r,f[i][j][oz+1].r,f[i][j][oz+2].r);
        df[i][j][oz].i = intplderiv(z[oz], 
           z[oz],z[oz+1],z[oz+2],f[i][j][oz].i,f[i][j][oz+1].i,f[i][j][oz+2].i);
        df[i][j][nz].r = intplderiv(z[nz], 
           z[nz],z[Nm1],z[Nm2],f[i][j][nz].r,f[i][j][Nm1].r,f[i][j][Nm2].r);
        df[i][j][nz].i = intplderiv(z[nz],
           z[nz],z[Nm1],z[Nm2],f[i][j][nz].i,f[i][j][Nm1].i,f[i][j][Nm2].i);
     }

     for (k=oz+1;k<nz;k++) {
        km1=k-1;  kp1=k+1;
        for (i=ox;i<=nx;i++)
        for (j=oy;j<=ny;j++) {   
           df[i][j][k].r = intplderiv(z[k], 
              z[km1],z[k],z[kp1],f[i][j][km1].r,f[i][j][k].r,f[i][j][kp1].r);
           df[i][j][k].i = intplderiv(z[k], 
              z[km1],z[k],z[kp1],f[i][j][km1].i,f[i][j][k].i,f[i][j][kp1].i);
        }
     }
  }
}


/* ============================================================================
   CFPDIFF3D:  Sept 6/93

   Find partial derivative of a complex 3D matrix.
   The derivative is taken of the index'th component of a spectral field.
     i.e.  diff( sum(n=0..nnx, A_n exp(I n alpha x)) ) =
              sum(n=0..nnx, I n alpha A_n exp(I n alpha x))

   ============================================================================
*/
void CFpdiff3D(f,ox,nx,oy,ny,oz,nz, kxyz,index, df)
   int ox,nx,oy,ny,oz,nz,index;
   double *kxyz;
   dcomplex ***f,***df;
{
  int n,m,l;

  if (index<1 || index>3) myerror("CFPDIFF3D: index out of bounds");

  if (index==1) {    /* take spectral derivative wrt x */
     for (n=ox;n<=nx;n++) 
     for (m=oy;m<=ny;m++) 
     for (l=oz;l<=nz;l++){
        df[n][m][l].r = -kxyz[n]*f[n][m][l].i;
        df[n][m][l].i = kxyz[n]*f[n][m][l].r;
     }
  } else if (index==2) {
     for (n=ox;n<=nx;n++)
     for (m=oy;m<=ny;m++) 
     for (l=oz;l<=nz;l++){
        df[n][m][l].r = -kxyz[m]*f[n][m][l].i;
        df[n][m][l].i = kxyz[m]*f[n][m][l].r;
     }
  } else if (index==3) {
     for (n=ox;n<=nx;n++)
     for (m=oy;m<=ny;m++) 
     for (l=oz;l<=nz;l++){
        df[n][m][l].r = -kxyz[l]*f[n][m][l].i;
        df[n][m][l].i = kxyz[l]*f[n][m][l].r;
     }
  }
}


/* ============================================================================
   ZzPDIFF3DO1:  Dec 22/99

   Find partial derivative with respect to index of a complex field defined 
   at evenly spaced grid points on a 3D matrix using a _first_ order finite 
   difference scheme.  Derivative is presumed to be calculated at the 
   intermediate gridpoint.   The direction of differentiation depends on
   the sign of index (i.e. see comments to PDIFF2DO1).
   By default, zeros are inserted in df where it is not defined.
   ============================================================================
*/
void Zzpdiff3DO1(f,ox,nx,oy,ny,oz,nz, dxyz,index, df)
   int ox,nx,oy,ny,oz,nz,index;
   double dxyz;
   zomplex ***f,***df;
{
  int aind;
  int i,j,k;
  int ip1,jp1,kp1,im1,jm1,km1;
  double idx,idy,idz;

  aind = iabs(index);
  if (aind==0 || aind>3) myerror("ZPDIFF3DO1: index out of bounds");

  if (index==1 && ox!=1) myerror("ZPDIFF3DO1: ox should equal one");
  if (index==2 && oy!=1) myerror("ZPDIFF3DO1: oy should equal one");
  if (index==3 && oz!=1) myerror("ZPDIFF3DO1: oz should equal one");
  if (index==-1 && ox!=0) myerror("ZPDIFF3DO1: ox should equal zero");
  if (index==-2 && oy!=0) myerror("ZPDIFF3DO1: oy should equal zero");
  if (index==-3 && oz!=0) myerror("ZPDIFF3DO1: oz should equal zero");

  if (index==1) {
     idx = 1.0/dxyz;
     for (i=1;i<nx;i++) {
        ip1=i+1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].re = idx*(f[ip1][j][k].re-f[i][j][k].re); 
           df[i][j][k].im = idx*(f[ip1][j][k].im-f[i][j][k].im); 
        }
     }
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) {
        df[0][j][k].re = df[1][j][k].re;
        df[0][j][k].im = df[1][j][k].im;
        df[nx][j][k].re = df[nx-1][j][k].re;
        df[nx][j][k].im = df[nx-1][j][k].im;
     }
  } else if (index==2) {
     idy = 1.0/dxyz;
     for (j=1;j<ny;j++) {
        jp1=j+1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].re = idy*(f[i][jp1][k].re-f[i][j][k].re); 
           df[i][j][k].im = idy*(f[i][jp1][k].im-f[i][j][k].im); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) {
        df[i][0][k].re = df[i][1][k].re;
        df[i][0][k].im = df[i][1][k].im;
        df[i][ny][k].re = df[i][ny-1][k].re;
        df[i][ny][k].im = df[i][ny-1][k].im;
     }
  } else if (index==3) {
     idz = 1.0/dxyz;
     for (k=1;k<nz;k++) {
        kp1=k+1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k].re = idz*(f[i][j][kp1].re-f[i][j][k].re); 
           df[i][j][k].im = idz*(f[i][j][kp1].im-f[i][j][k].im); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (j=oy;j<=ny;j++) {
        df[i][j][0].re = df[i][j][1].re;
        df[i][j][0].im = df[i][j][1].im;
        df[i][j][nz].re = df[i][j][nz-1].re;
        df[i][j][nz].im = df[i][j][nz-1].im;
     }
  } else if (index==-1) {
     idx = 1.0/dxyz;
     for (i=1;i<=nx;i++) {
        im1=i-1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].re = idx*(f[i][j][k].re-f[im1][j][k].re); 
           df[i][j][k].im = idx*(f[i][j][k].im-f[im1][j][k].im); 
        }
     }
  } else if (index==-2) {
     idy = 1.0/dxyz;
     for (j=1;j<=ny;j++) {
        jm1=j-1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].re = idy*(f[i][j][k].re-f[i][jm1][k].re); 
           df[i][j][k].im = idy*(f[i][j][k].im-f[i][jm1][k].im); 
        }
     }
  } else if (index==-3) {
     idz = 1.0/dxyz;
     for (k=1;k<=nz;k++) {
        km1=k-1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k].re = idz*(f[i][j][k].re-f[i][j][km1].re); 
           df[i][j][k].im = idz*(f[i][j][k].im-f[i][j][km1].im); 
        }
     }
  }
}

/* ============================================================================
   ZPPDIFF3DO1:  Dec 7/93

   Find partial derivative with respect to index of a complex field defined 
   at evenly spaced grid points on a 3D matrix using a _first_ order finite 
   difference scheme.  Derivative is presumed to be calculated at the 
   intermediate gridpoint.   The direction of differentiation depends on
   the sign of index (i.e. see comments to PDIFF2DO1).
   By default, zeros are inserted in df where it is not defined.

   This procedure differs from ZPDIFFO1 in that is assumes the domain
   has periodic boundaries.
   ============================================================================
*/

void ZPpdiff3DO1(f,ox,nx,oy,ny,oz,nz, dxyz,index, df)
   int ox,nx,oy,ny,oz,nz,index;
   double dxyz;
   zomplex ***f,***df;
{
  int aind;
  int i,j,k;
  int ip1,jp1,kp1,im1,jm1,km1;
  double idx,idy,idz;

  aind = iabs(index);
  if (aind==0 || aind>3) myerror("ZPPDIFF3DO1: index out of bounds");

  if (aind==1 && ox!=0) myerror("ZPPDIFF3DO1: ox should equal zero");
  if (aind==2 && oy!=0) myerror("ZPPDIFF3DO1: oy should equal zero");
  if (aind==3 && oz!=0) myerror("ZPPDIFF3DO1: oz should equal zero");

  if (index==1) {
     idx = 1.0/dxyz;
     for (i=0;i<nx;i++) {
        ip1=i+1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].re = idx*(f[ip1][j][k].re-f[i][j][k].re); 
           df[i][j][k].im = idx*(f[ip1][j][k].im-f[i][j][k].im); 
        }
     }
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) df[nx][j][k] = df[0][j][k];
  } else if (index==2) {
     idy = 1.0/dxyz;
     for (j=0;j<ny;j++) {
        jp1=j+1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].re = idy*(f[i][jp1][k].re-f[i][j][k].re); 
           df[i][j][k].im = idy*(f[i][jp1][k].im-f[i][j][k].im); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) df[i][ny][k] = df[i][0][k];
  } else if (index==3) {
     idz = 1.0/dxyz;
     for (k=0;k<nz;k++) {
        kp1=k+1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k].re = idz*(f[i][j][kp1].re-f[i][j][k].re); 
           df[i][j][k].im = idz*(f[i][j][kp1].im-f[i][j][k].im); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (j=oy;j<=ny;j++) df[i][j][nz] = df[i][j][0];
  } else if (index==-1) {
     idx = 1.0/dxyz;
     for (i=1;i<=nx;i++) {
        im1=i-1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].re = idx*(f[i][j][k].re-f[im1][j][k].re); 
           df[i][j][k].im = idx*(f[i][j][k].im-f[im1][j][k].im); 
        }
     }
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) df[0][j][k] = df[nx][j][k];
  } else if (index==-2) {
     idy = 1.0/dxyz;
     for (j=1;j<=ny;j++) {
        jm1=j-1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].re = idy*(f[i][j][k].re-f[i][jm1][k].re); 
           df[i][j][k].im = idy*(f[i][j][k].im-f[i][jm1][k].im); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) df[i][0][k] = df[i][ny][k];
  } else if (index==-3) {
     idz = 1.0/dxyz;
     for (k=1;k<=nz;k++) {
        km1=k-1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k].re = idz*(f[i][j][k].re-f[i][j][km1].re); 
           df[i][j][k].im = idz*(f[i][j][k].im-f[i][j][km1].im); 
        }
     }
     for (i=ox;i<=nx;i++)  
     for (j=oy;j<=ny;j++) df[i][j][0] = df[i][j][nz];
  }
}

/* ============================================================================
   ZPDIFF3D:  Sept 20/93

   Find partial derivative of a complex 3D matrix with respect to index
   using a second order finite difference scheme.  Function is defined
   on an evenly spaced grid.
   ============================================================================
*/

void Zpdiff3D(f,ox,nx,oy,ny,oz,nz, dxyz,index, df)
   int ox,nx,oy,ny,oz,nz,index;
   double dxyz;
   zomplex ***f,***df;
{
  int i,j,k;
  int ip1,im1,jp1,jm1,kp1,km1,Nm1,Nm2;
  double idx2,idy2,idz2;

  if (index<1 || index>3) myerror("ZPDIFF3D: index out of bounds");

  if (index==1) {
     idx2 = 0.5/dxyz;
     Nm1=nx-1;   Nm2=nx-2;
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) { 
        df[ox][j][k].re = 
           idx2*(-3*f[ox][j][k].re + 4*f[ox+1][j][k].re - f[ox+2][j][k].re);
        df[ox][j][k].im = 
           idx2*(-3*f[ox][j][k].im + 4*f[ox+1][j][k].im - f[ox+2][j][k].im);
        df[nx][j][k].re = 
           idx2*(3*f[nx][j][k].re - 4*f[Nm1][j][k].re + f[Nm2][j][k].re);
        df[nx][j][k].im = 
           idx2*(3*f[nx][j][k].im - 4*f[Nm1][j][k].im + f[Nm2][j][k].im);
     }
     for (i=ox+1;i<nx;i++) {
        ip1=i+1;  im1=i-1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].re = idx2*(f[ip1][j][k].re-f[im1][j][k].re); 
           df[i][j][k].im = idx2*(f[ip1][j][k].im-f[im1][j][k].im); 
        }
     }
  } else if (index==2) {
     idy2 = 0.5/dxyz;
     Nm1=ny-1;   Nm2=ny-2;
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) { 
        df[i][oy][k].re = 
           idy2*(-3*f[i][oy][k].re + 4*f[i][oy+1][k].re - f[i][oy+2][k].re);
        df[i][oy][k].im = 
           idy2*(-3*f[i][oy][k].im + 4*f[i][oy+1][k].im - f[i][oy+2][k].im);
        df[i][ny][k].re = 
           idy2*(3*f[i][ny][k].re - 4*f[i][Nm1][k].re + f[i][Nm2][k].re);
        df[i][ny][k].im = 
           idy2*(3*f[i][ny][k].im - 4*f[i][Nm1][k].im + f[i][Nm2][k].im);
     }
     for (j=oy+1;j<ny;j++) {
        jp1=j+1;  jm1=j-1;
        for (i=ox;i<=nx;i++)  
        for (k=ox;k<=nz;k++) { 
           df[i][j][k].re = idy2*(f[i][jp1][k].re-f[i][jm1][k].re); 
           df[i][j][k].im = idy2*(f[i][jp1][k].im-f[i][jm1][k].im); 
        }
     }
  } else if (index==3) {
     idz2 = 0.5/dxyz;
     Nm1=nz-1;   Nm2=nz-2;
     for (i=ox;i<=nx;i++)  
     for (j=oy;j<=ny;j++) {  
        df[i][j][oz].re = 
           idz2*(-3*f[i][j][oz].re + 4*f[i][j][oz+1].re - f[i][j][oz+2].re);
        df[i][j][oz].im = 
           idz2*(-3*f[i][j][oz].im + 4*f[i][j][oz+1].im - f[i][j][oz+2].im);
        df[i][j][nz].re = 
           idz2*(3*f[i][j][nz].re - 4*f[i][j][Nm1].re + f[i][j][Nm2].re);
        df[i][j][nz].im = 
           idz2*(3*f[i][j][nz].im - 4*f[i][j][Nm1].im + f[i][j][Nm2].im);
     }
     for (k=oz+1;k<nz;k++) {
        kp1=k+1; km1=k-1;
        for (i=ox;i<=nx;i++)  
        for (j=oy;j<=ny;j++) {  
           df[i][j][k].re = idz2*(f[i][j][kp1].re-f[i][j][km1].re); 
           df[i][j][k].im = idz2*(f[i][j][kp1].im-f[i][j][km1].im); 
        }
     }
  }
}

/* ============================================================================
   ZPDIFF3D_gen:  Sept 20/93

   Find partial derivative of a complex 3D matrix with respect to index
   using a second order finite difference scheme.  Function is not
   necessarily defined at evenly space points.
   ============================================================================
*/

void Zpdiff3D_gen(f, x,ox,nx, y,oy,ny, z,oz,nz, index, df)
   int ox,nx,oy,ny,oz,nz,index;
   double *x,*y,*z;
   zomplex ***f,***df;
{
  int i,j,k;
  int ip1,im1,jp1,jm1,kp1,km1,Nm1,Nm2;

  if (index<1 || index>3) myerror("ZPDIFF3D_gen: index out of bounds");

  if (index==1) {
     Nm1=nx-1; Nm2=nx-2;
     for (j=oy;j<=ny;j++)  
     for (k=oz;k<=nz;k++) { 
        df[ox][j][k].re = intplderiv(x[ox], 
           x[ox],x[ox+1],x[ox+2],f[ox][j][k].re,f[ox+1][j][k].re,f[ox+2][j][k].re);
        df[ox][j][k].im = intplderiv(x[ox], 
           x[ox],x[ox+1],x[ox+2],f[ox][j][k].im,f[ox+1][j][k].im,f[ox+2][j][k].im);
        df[nx][j][k].re = intplderiv(x[nx], 
           x[nx],x[Nm1],x[Nm2],f[nx][j][k].re,f[Nm1][j][k].re,f[Nm2][j][k].re);
        df[nx][j][k].im = intplderiv(x[nx],
           x[nx],x[Nm1],x[Nm2],f[nx][j][k].im,f[Nm1][j][k].im,f[Nm2][j][k].im);
     }

     for (i=ox+1;i<nx;i++) {  
        im1=i-1; ip1=i+1;
        for (j=oy;j<=ny;j++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].re = intplderiv(x[i], 
              x[im1],x[i],x[ip1],f[im1][j][k].re,f[i][j][k].re,f[ip1][j][k].re);
           df[i][j][k].im = intplderiv(x[i], 
              x[im1],x[i],x[ip1],f[im1][j][k].im,f[i][j][k].im,f[ip1][j][k].im);
        }
     }

  } else if (index==2) {
     Nm1=ny-1; Nm2=ny-2;
     for (i=ox;i<=nx;i++)  
     for (k=oz;k<=nz;k++) { 
        df[i][oy][k].re = intplderiv(y[oy], 
           y[oy],y[oy+1],y[oy+2],f[i][oy][k].re,f[i][oy+1][k].re,f[i][oy+2][k].re);
        df[i][oy][k].im = intplderiv(y[oy], 
           y[oy],y[oy+1],y[oy+2],f[i][oy][k].im,f[i][oy+1][k].im,f[i][oy+2][k].im);
        df[i][ny][k].re = intplderiv(y[ny], 
           y[ny],y[Nm1],y[Nm2],f[i][ny][k].re,f[i][Nm1][k].re,f[i][Nm2][k].re);
        df[i][ny][k].im = intplderiv(y[ny],
           y[ny],y[Nm1],y[Nm2],f[i][ny][k].im,f[i][Nm1][k].im,f[i][Nm2][k].im);
     }

     for (j=oy+1;j<ny;j++) {  
        jm1=j-1;  jp1=j+1;
        for (i=ox;i<=nx;i++)  
        for (k=oz;k<=nz;k++) { 
           df[i][j][k].re = intplderiv(y[j], 
              y[jm1],y[j],y[jp1],f[i][jm1][k].re,f[i][j][k].re,f[i][jp1][k].re);
           df[i][j][k].im = intplderiv(y[j], 
              y[jm1],y[j],y[jp1],f[i][jm1][k].im,f[i][j][k].im,f[i][jp1][k].im);
        }
     }

  } else if (index==3) {
     Nm1=nz-1;  Nm2=nz-2;
     for (i=ox;i<=nx;i++)
     for (j=oy;j<=ny;j++) {   
        df[i][j][oz].re = intplderiv(z[oz], 
           z[oz],z[oz+1],z[oz+2],f[i][j][oz].re,f[i][j][oz+1].re,f[i][j][oz+2].re);
        df[i][j][oz].im = intplderiv(z[oz], 
           z[oz],z[oz+1],z[oz+2],f[i][j][oz].im,f[i][j][oz+1].im,f[i][j][oz+2].im);
        df[i][j][nz].re = intplderiv(z[nz], 
           z[nz],z[Nm1],z[Nm2],f[i][j][nz].re,f[i][j][Nm1].re,f[i][j][Nm2].re);
        df[i][j][nz].im = intplderiv(z[nz],
           z[nz],z[Nm1],z[Nm2],f[i][j][nz].im,f[i][j][Nm1].im,f[i][j][Nm2].im);
     }

     for (k=oz+1;k<nz;k++) {
        km1=k-1;  kp1=k+1;
        for (i=ox;i<=nx;i++)
        for (j=oy;j<=ny;j++) {   
           df[i][j][k].re = intplderiv(z[k], 
              z[km1],z[k],z[kp1],f[i][j][km1].re,f[i][j][k].re,f[i][j][kp1].re);
           df[i][j][k].im = intplderiv(z[k], 
              z[km1],z[k],z[kp1],f[i][j][km1].im,f[i][j][k].im,f[i][j][kp1].im);
        }
     }
  }
}


/* ============================================================================
   ZFPDIFF3D:  Sept 6/93

   Find partial derivative of a complex 3D matrix.
   The derivative is taken of the index'th component of a spectral field.
     i.e.  diff( sum(n=0..nnx, A_n exp(I n alpha x)) ) =
              sum(n=0..nnx, I n alpha A_n exp(I n alpha x))

   ============================================================================
*/
void ZFpdiff3D(f,ox,nx,oy,ny,oz,nz, kxyz,index, df)
   int ox,nx,oy,ny,oz,nz,index;
   double *kxyz;
   zomplex ***f,***df;
{
  int n,m,l;

  if (index<1 || index>3) myerror("ZFPDIFF3D: index out of bounds");

  if (index==1) {    /* take spectral derivative wrt x */
     for (n=ox;n<=nx;n++) 
     for (m=oy;m<=ny;m++) 
     for (l=oz;l<=nz;l++){
        df[n][m][l].re = -kxyz[n]*f[n][m][l].im;
        df[n][m][l].im = kxyz[n]*f[n][m][l].re;
     }
  } else if (index==2) {
     for (n=ox;n<=nx;n++)
     for (m=oy;m<=ny;m++) 
     for (l=oz;l<=nz;l++){
        df[n][m][l].re = -kxyz[m]*f[n][m][l].im;
        df[n][m][l].im = kxyz[m]*f[n][m][l].re;
     }
  } else if (index==3) {
     for (n=ox;n<=nx;n++)
     for (m=oy;m<=ny;m++) 
     for (l=oz;l<=nz;l++){
        df[n][m][l].re = -kxyz[l]*f[n][m][l].im;
        df[n][m][l].im = kxyz[l]*f[n][m][l].re;
     }
  }
}


