/* =======================================================================
   Various routines to perform real and complex arithmetic on 
   2D and 3D matrices.


   Functions:

   Cadd2D     _Cadd2D
   Csub2D     _Csub2D
   Cmul2D     
   Ccmul2D    _Ccmul2D
   CMmul2D     

   Radd3D     _Radd3D
   Rsub3D     _Rsub3D
   Rmul3D     _Rmul3D
   Rdiv3D     _Rdiv3D
   Rcmul3D    _Rcmul3D
   Rmul3D2D   _Rmul3D2D
   Cadd3D     _Cadd3D
   Csub3D     _Csub3D
   Ccmul3D    _Ccmul3D

   CRvmul3D   _CRvmul3D

   MBsmooth
   MBzero
   =======================================================================
*/
#include "complex.h"
#include "alloc_space.h"

void _Ccmul2D(f,c, ox,nx,oy,ny)
   int ox,nx,oy,ny;
   double c;
   dcomplex **f;
{
   int i,j;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++) {
      f[i][j].r *= c;
      f[i][j].i *= c;
   }
}

void _Cadd2D(f1,f2, ox,nx,oy,ny)
   int ox,nx,oy,ny;
   dcomplex **f1,**f2;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++) {
      f1[i][j].r += f2[i][j].r;
      f1[i][j].i += f2[i][j].i;
   }
}

void _Csub2D(f1,f2, ox,nx,oy,ny)
   int ox,nx,oy,ny;
   dcomplex **f1,**f2;
{
   int i,j;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++) {
      f1[i][j].r -= f2[i][j].r;
      f1[i][j].i -= f2[i][j].i;
   }
}

void Radd3D(f1,f2, ox,nx,oy,ny,oz,nz, f3)
   int ox,nx,oy,ny,oz,nz;
   double ***f1,***f2,***f3;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++)
   for (k=oz;k<=nz;k++) f3[i][j][k] = f1[i][j][k]+f2[i][j][k];
}

void Cmul2D(f1,f2, ox,nx,oy,ny, f3)
   int ox,nx,oy,ny;
   dcomplex **f1,**f2,**f3;
{
   int i,j;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++) f3[i][j] = Cmul(f1[i][j],f2[i][j]);
}

void CMmul2D(f1,f2, ox,nx,oy,ny,oz,nz, f3)
   int ox,nx,oy,ny,oz,nz;
   dcomplex **f1,**f2,**f3;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (k=oz;k<=nz;k++) {
      f3[i][k].r = f3[i][k].i = 0.0;
      for (j=oy;j<=ny;j++) { 
         f3[i][k].r += f1[i][j].r*f2[j][k].r - f1[i][j].i*f2[i][j].i;
         f3[i][k].i += f1[i][j].r*f2[j][k].i + f1[i][j].i*f2[i][j].r;
      }
   }
}

void _Radd3D(f1,f2, ox,nx,oy,ny,oz,nz)
   int ox,nx,oy,ny,oz,nz;
   double ***f1,***f2;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++)
   for (k=oz;k<=nz;k++) f1[i][j][k] += f2[i][j][k];
}

void Rsub3D(f1,f2, ox,nx,oy,ny,oz,nz, f3)
   int ox,nx,oy,ny,oz,nz;
   double ***f1,***f2,***f3;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++)
   for (k=oz;k<=nz;k++) f3[i][j][k] = f1[i][j][k]-f2[i][j][k];
}

void _Rsub3D(f1,f2, ox,nx,oy,ny,oz,nz)
   int ox,nx,oy,ny,oz,nz;
   double ***f1,***f2;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++)
   for (k=oz;k<=nz;k++) f1[i][j][k] -= f2[i][j][k];
}

void Rmul3D(f1,f2, ox,nx,oy,ny,oz,nz, f3)
   int ox,nx,oy,ny,oz,nz;
   double ***f1,***f2,***f3;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++)
   for (k=oz;k<=nz;k++) f3[i][j][k] = f1[i][j][k]*f2[i][j][k];
}

void _Rmul3D(f1,f2, ox,nx,oy,ny,oz,nz)
   int ox,nx,oy,ny,oz,nz;
   double ***f1,***f2;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++)
   for (k=oz;k<=nz;k++) f1[i][j][k] *= f2[i][j][k];
}

void Rdiv3D(f1,f2, ox,nx,oy,ny,oz,nz, f3)
   int ox,nx,oy,ny,oz,nz;
   double ***f1,***f2,***f3;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++)
   for (k=oz;k<=nz;k++) f3[i][j][k] = f1[i][j][k]/f2[i][j][k];
}

void _Rdiv3D(f1,f2, ox,nx,oy,ny,oz,nz)
   int ox,nx,oy,ny,oz,nz;
   double ***f1,***f2;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++)
   for (k=oz;k<=nz;k++) f1[i][j][k] /= f2[i][j][k];
}

void Rcmul3D(f1,c, ox,nx,oy,ny,oz,nz, f2)
   int ox,nx,oy,ny,oz,nz;
   double c;
   double ***f1,***f2;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++)
   for (k=oz;k<=nz;k++) f2[i][j][k] = c*f1[i][j][k];
}

void _Rcmul3D(f1,c, ox,nx,oy,ny,oz,nz)
   int ox,nx,oy,ny,oz,nz;
   double c;
   double ***f1;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++)
   for (k=oz;k<=nz;k++) f1[i][j][k] *= c;
}

void Rmul3D1D(f1,f2,index, ox,nx,oy,ny,oz,nz, f3)
   int index,ox,nx,oy,ny,oz,nz;
   double *f2;
   double ***f1,***f3;
{
   int i,j,k;

   if (index==1)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f3[i][j][k] = f1[i][j][k]*f2[i];
   else if (index==2)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f3[i][j][k] = f1[i][j][k]*f2[j];
   else if (index==3)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f3[i][j][k] = f1[i][j][k]*f2[k];
   else 
      myerror("RMUL3D1D: index out of bounds");
}

void _Rmul3D1D(f1,f2,index, ox,nx,oy,ny,oz,nz)
   int index,ox,nx,oy,ny,oz,nz;
   double *f2;
   double ***f1;
{
   int i,j,k;

   if (index==1)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f1[i][j][k] *= f2[i];
   else if (index==2)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f1[i][j][k] *= f2[j];
   else if (index==3)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f1[i][j][k] *= f2[k];
   else 
      myerror("_RMUL3D1D: index out of bounds");
}

void Rmul3D2D(f1,f2,index, ox,nx,oy,ny,oz,nz, f3)
   int index,ox,nx,oy,ny,oz,nz;
   double **f2;
   double ***f1,***f3;
{
   int i,j,k;

   if (index==1)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f3[i][j][k] = f1[i][j][k]*f2[j][k];
   else if (index==2)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f3[i][j][k] = f1[i][j][k]*f2[i][k];
   else if (index==3)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f3[i][j][k] = f1[i][j][k]*f2[i][j];
   else 
      myerror("RMUL3D2D: index out of bounds");
}

void _Rmul3D2D(f1,f2,index, ox,nx,oy,ny,oz,nz)
   int index,ox,nx,oy,ny,oz,nz;
   double **f2;
   double ***f1;
{
   int i,j,k;

   if (index==1)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f1[i][j][k] *= f2[j][k];
   else if (index==2)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f1[i][j][k] *= f2[i][k];
   else if (index==3)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f1[i][j][k] *= f2[i][j];
   else 
      myerror("_RMUL3D2D: index out of bounds");
}

void Rdiv3D2D(f1,f2,index, ox,nx,oy,ny,oz,nz, f3)
   int index,ox,nx,oy,ny,oz,nz;
   double **f2;
   double ***f1,***f3;
{
   int i,j,k;

   if (index==1)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f3[i][j][k] = f1[i][j][k]/f2[j][k];
   else if (index==2)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f3[i][j][k] = f1[i][j][k]/f2[i][k];
   else if (index==3)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f3[i][j][k] = f1[i][j][k]/f2[i][j];
   else 
      myerror("RDIV3D2D: index out of bounds");
}

void _Rdiv3D2D(f1,f2,index, ox,nx,oy,ny,oz,nz)
   int index,ox,nx,oy,ny,oz,nz;
   double **f2;
   double ***f1;
{
   int i,j,k;

   if (index==1)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f1[i][j][k] /= f2[j][k];
   else if (index==2)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f1[i][j][k] /= f2[i][k];
   else if (index==3)
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) f1[i][j][k] /= f2[i][j];
   else 
      myerror("_RDIV3D2D: index out of bounds");
}

void _Ccmul3D(f,c, ox,nx,oy,ny,oz,nz)
   int ox,nx,oy,ny,oz,nz;
   double c;
   dcomplex ***f;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++)
   for (k=oz;k<=nz;k++) { 
      f[i][j][k].r *= c;
      f[i][j][k].i *= c;
   }
}

void _CRvmul3D(f,v,index, ox,nx,oy,ny,oz,nz)
   int index,ox,nx,oy,ny,oz,nz;
   double *v;
   dcomplex ***f;
{
   int i,j,k;

   if (index==1) {
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) { 
         f[i][j][k].r *= v[i];
         f[i][j][k].i *= v[i];
      }
   } else if (index==2) {
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) { 
         f[i][j][k].r *= v[j];
         f[i][j][k].i *= v[j];
      }
   } else if (index==3) {
      for (i=ox;i<=nx;i++) 
      for (j=oy;j<=ny;j++)
      for (k=oz;k<=nz;k++) { 
         f[i][j][k].r *= v[k];
         f[i][j][k].i *= v[k];
      }
   }
}

void _Cadd3D(f1,f2, ox,nx,oy,ny,oz,nz)
   int ox,nx,oy,ny,oz,nz;
   dcomplex ***f1,***f2;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++)
   for (k=oz;k<=nz;k++) { 
      f1[i][j][k].r += f2[i][j][k].r;
      f1[i][j][k].i += f2[i][j][k].i;
   }
}

void _Csub3D(f1,f2, ox,nx,oy,ny,oz,nz)
   int ox,nx,oy,ny,oz,nz;
   dcomplex ***f1,***f2;
{
   int i,j,k;

   for (i=ox;i<=nx;i++) 
   for (j=oy;j<=ny;j++)
   for (k=oz;k<=nz;k++) { 
      f1[i][j][k].r -= f2[i][j][k].r;
      f1[i][j][k].i -= f2[i][j][k].i;
   }
}

/* =========================================================================
   RMBSMOOTH2D: FEB 15/94
 
   Take a real field and smooth out the meridional boundaries by
   linear extrapolation.  This routine is used following a the calculation
   of the diffusion and heat dissipation contribution to the time
   derivatives.
   =========================================================================
*/
void Rmbsmooth2D(f, ox,nx,oy,ny, sf)
   int ox,nx,oy,ny;
   double **f,**sf;
{
  int i,j;

  for (i=ox;i<=nx;i++) 
  for (j=oy+1;j<ny;j++) sf[i][j] = f[i][j];

  for (i=ox;i<=nx;i++) {
     sf[i][oy] = 2*f[i][oy+1] - f[i][oy+2];
     sf[i][ny] = 2*f[i][ny-1] - f[i][ny-2];
  }
}

/* =========================================================================
   _RMBSMOOTH2D: FEB 15/94
 
   Take a real field and smooth out the meridional boundaries by
   linear extrapolation.  This routine is used following a the calculation
   of the diffusion and heat dissipation contribution to the time
   derivatives.
   =========================================================================
*/
void _Rmbsmooth2D(f, ox,nx,oy,ny)
   int ox,nx,oy,ny;
   double **f;
{
  int n,k;

  for (n=ox;n<=nx;n++) {
     f[n][oy] = 2*f[n][oy+1] - f[n][oy+2];
     f[n][ny] = 2*f[n][ny-1] - f[n][ny-2];
  }
}

/* =========================================================================
   RMBSMOOTH: FEB 15/94
 
   Take a real field and smooth out the meridional boundaries by
   linear extrapolation.  This routine is used following a the calculation
   of the diffusion and heat dissipation contribution to the time
   derivatives.
   =========================================================================
*/
void Rmbsmooth(f, ox,nx,oy,ny,os,ns, sf)
   int ox,nx,oy,ny,os,ns;
   double ***f,***sf;
{
  int i,j,k;

  for (i=ox;i<=nx;i++)
  for (j=oy+1;j<ny;j++)
  for (k=os;k<=ns;k++) sf[i][j][k] = f[i][j][k];

  for (i=ox;i<=nx;i++)
  for (k=os;k<=ns;k++) {
     sf[i][oy][k] = 2*f[i][oy+1][k] - f[i][oy+2][k];
     sf[i][ny][k] = 2*f[i][ny-1][k] - f[i][ny-2][k];
  }
}

/* =========================================================================
   _RMBSMOOTH: FEB 15/94
 
   Take a real field and smooth out the meridional boundaries by
   linear extrapolation.  This routine is used following a the calculation
   of the diffusion and heat dissipation contribution to the time
   derivatives.
   =========================================================================
*/
void _Rmbsmooth(f, ox,nx,oy,ny,os,ns)
   int ox,nx,oy,ny,os,ns;
   double ***f;
{
  int n,k;

  for (n=ox;n<=nx;n++)
  for (k=os;k<=ns;k++) {
     f[n][oy][k] = 2*f[n][oy+1][k] - f[n][oy+2][k];
     f[n][ny][k] = 2*f[n][ny-1][k] - f[n][ny-2][k];
  }
}

/* =========================================================================
   _CMBSMOOTH2D: NOV 25/93
 
   Take a complex field and smooth out the meridional boundaries by
   linear extrapolation.  This routine is used following a the calculation
   of the diffusion and heat dissipation contribution to the time
   derivatives.
   =========================================================================
*/
void _Cmbsmooth2D(f, ox,nnx,oy,ny)
   int ox,nnx,oy,ny;
   dcomplex **f;
{
  int n;

  for (n=ox;n<=nnx;n++) {
     f[n][oy].r = 2*f[n][oy+1].r - f[n][oy+2].r;
     f[n][oy].i = 2*f[n][oy+1].i - f[n][oy+2].i;
     f[n][ny].r = 2*f[n][ny-1].r - f[n][ny-2].r;
     f[n][ny].i = 2*f[n][ny-1].i - f[n][ny-2].i;
  }
}

/* =========================================================================
   _CMBSMOOTH2D_2: MAY 3/94
 
   Take a complex field and smooth out the meridional boundaries by
   linear extrapolation.  This routine differs from _CMBSMOOTH2D in
   that smoothing takes place over _two_ layers near boundaries.
   =========================================================================
*/
void _Cmbsmooth2D_2(f, ox,nnx,oy,ny)
   int ox,nnx,oy,ny;
   dcomplex **f;
{
  int n;

  for (n=ox;n<=nnx;n++) {
     f[n][oy+1].r = 2*f[n][oy+2].r - f[n][oy+3].r;
     f[n][oy+1].i = 2*f[n][oy+2].i - f[n][oy+3].i;
     f[n][oy].r = 2*f[n][oy+1].r - f[n][oy+2].r;
     f[n][oy].i = 2*f[n][oy+1].i - f[n][oy+2].i;
     f[n][ny-1].r = 2*f[n][ny-2].r - f[n][ny-3].r;
     f[n][ny-1].i = 2*f[n][ny-2].i - f[n][ny-3].i;
     f[n][ny].r = 2*f[n][ny-1].r - f[n][ny-2].r;
     f[n][ny].i = 2*f[n][ny-1].i - f[n][ny-2].i;
  }
}

/* =========================================================================
   _CMBSMOOTH: NOV 20/93
 
   Take a complex field and smooth out the meridional boundaries by
   linear extrapolation.  This routine is used following a the calculation
   of the diffusion and heat dissipation contribution to the time
   derivatives.
   =========================================================================
*/
void _Cmbsmooth(f, ox,nnx,oy,ny,os,ns)
   int ox,nnx,oy,ny,os,ns;
   dcomplex ***f;
{
  int n,k;

  for (n=ox;n<=nnx;n++)
  for (k=os;k<=ns;k++) {
     f[n][oy][k].r = 2*f[n][oy+1][k].r - f[n][oy+2][k].r;
     f[n][oy][k].i = 2*f[n][oy+1][k].i - f[n][oy+2][k].i;
     f[n][ny][k].r = 2*f[n][ny-1][k].r - f[n][ny-2][k].r;
     f[n][ny][k].i = 2*f[n][ny-1][k].i - f[n][ny-2][k].i;
  }
}

/* =========================================================================
   _RMBZERO2D: Dec 15/93
 
   Take a 2D real field and zero the meridional boundaries (except DC mode).
   This routine is applied to the v and time derivative fields.
   =========================================================================
*/
void _Rmbzero2D(f, ox,nx,oy,ny)
   int ox,nx,oy,ny;
   double **f;
{
  int n;

  for (n=ox;n<=nx;n++) f[n][oy] = f[n][ny] = 0.0;
}

/* =========================================================================
   _RMBZERO: Dec 15/93
 
   Take a real field and zero the meridional boundaries (except DC mode).
   This routine is applied to the v and time derivative fields.
   =========================================================================
*/
void _Rmbzero(f, ox,nx,oy,ny,os,ns)
   int ox,nx,oy,ny,os,ns;
   double ***f;
{
  int n,k;

  for (n=ox;n<=nx;n++)
  for (k=os;k<=ns;k++) f[n][oy][k] = f[n][ny][k] = 0.0;
}

/* =========================================================================
   _CMBZERO2D: NOV 25/93
 
   Take a 2D complex field and zero the meridional boundaries (except DC mode).
   This routine is applied to the v and time derivative fields.
   =========================================================================
*/
void _Cmbzero2D(f, ox,nnx,oy,ny)
   int ox,nnx,oy,ny;
   dcomplex **f;
{
  int n;
  dcomplex cZero;

  cZero = Complex(0.0,0.0);

  for (n=ox;n<=nnx;n++) f[n][oy] = f[n][ny] = cZero;
}

/* =========================================================================
   _CMBZERO: NOV 20/93
 
   Take a complex field and zero the meridional boundaries (except DC mode).
   This routine is applied to the v and time derivative fields.
   =========================================================================
*/
void _Cmbzero(f, ox,nnx,oy,ny,os,ns)
   int ox,nnx,oy,ny,os,ns;
   dcomplex ***f;
{
  int n,k;
  dcomplex cZero;

  cZero = Complex(0.0,0.0);

  for (n=ox;n<=nnx;n++)
  for (k=os;k<=ns;k++) f[n][oy][k] = f[n][ny][k] = cZero;
}

