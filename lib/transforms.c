/* ===================================================================
   Routines to take transforms of special functions. 
   See myfft for (fast) Fourier Transforms.
   ===================================================================
 */
#include <stdio.h>
#include <math.h>
#include "complex.h"
#include "alloc_space.h"
#include "specialfn.h"
#include "macros.h"

/* -------------------------------------------------------------
   DEFJ0ARRAY: Feb 13, 2004
   
   Define matrix of tabulated values of J0(z_i r_j) 
     z_i is i'th zero (i=1..nz)
     r_i is radius lying between 0 and 1  (j=0 .. nr)
   -------------------------------------------------------------
*/
void defj0array(J0array,nr,nz,J0z)
   int nr,nz;
   double *J0z;
   double **J0array;
{
  int i,j;
  double dr,r;

  dr = 1.0/((double) nr);
  r = 0.0;
  for (j=0;j<=nr;j++) {
    for (i=1;i<=nz;i++) J0array[i][j] =  j0(J0z[i]*r);
    r += dr;
  }
}

/* -------------------------------------------------------------
   DEFJ1ARRAY: Feb 18, 2004
   
   Define matrix of tabulated values of J1(z_i r_j) 
     z_i is i'th zero (i=1..nz)
     r_i is radius lying between 0 and 1  (j=0 .. nr)
   -------------------------------------------------------------
*/
void defj1array(J1array,nr,nz,J1z)
   int nr,nz;
   double *J1z;
   double **J1array;
{
  int i,j;
  double dr,r;

  dr = 1.0/((double) nr);
  r = 0.0;
  for (j=0;j<=nr;j++) {
    for (i=1;i<=nz;i++) J1array[i][j] =  j1(J1z[i]*r);
    r += dr;
  }
}

/* -------------------------------------------------------------
   DEFJNARRAY: Jun 9, 2012
   
   Define matrix of tabulated values of J_m(z_i r_j) 
     z_i is i'th zero (i=1..nz)
     r_i is radius lying between 0 and 1  (j=0 .. nr)
   -------------------------------------------------------------
*/
void defjnarray(JNarray,n,nr,nz,JNz)
   int n,nr,nz;
   double *JNz;
   double **JNarray;
{
  int i,j;
  double dr,r;

  dr = 1.0/((double) nr);
  r = 0.0;
  for (j=0;j<=nr;j++) {
    for (i=1;i<=nz;i++) JNarray[i][j] =  jn(n,JNz[i]*r);
    r += dr;
  }
}

/* -------------------------------------------------------------
   NRMJ0SER: Feb 13, 2004
   
   Given zeros J0z of J0 Bessel function, define normalization of 
   J0 Bessel series: nrmJ0 = 1/([J1(J0z)]^2/2).
   (This is equivalent to 1/int_0^1 [J0(J0z r)]^2 r dr)
   -------------------------------------------------------------
*/
void nrmJ0ser(nrmJ0,nz,J0z)
   int nz;
   double *nrmJ0,*J0z;
{
  int i;
  double j1v;

  for (i=1;i<=nz;i++) {
    j1v = j1(J0z[i]);
    nrmJ0[i] = 2.0/sqr(j1v);
  }
}

/* -------------------------------------------------------------
   NRMJ1SER: Feb 18, 2004
   
   Given zeros J1z of J1 Bessel function, define normalization of 
   J1 Bessel series: nrmJ1 = 1/([J2(J1z)]^2/2).   [CHECK!!!]
   (This is equivalent to 1/int_0^1 [J1(J1z r)]^2 r dr)
   -------------------------------------------------------------
*/
void nrmJ1ser(nrmJ1,nz,J1z)
   int nz;
   double *nrmJ1,*J1z;
{
  int i;
  double j2v;

  for (i=1;i<=nz;i++) {
    j2v = jn(2,J1z[i]);
    nrmJ1[i] = 2.0/sqr(j2v);
  }
}

/* -------------------------------------------------------------
   NRMJNSER: Jun 9, 2012
   
   Given zeros Jnz of Jn Bessel function, define normalization of 
   Jn Bessel series: nrmJn = 1/([J(n+1)(Jnz)]^2/2).   [CHECK!!!]
   (This is equivalent to 1/int_0^1 [Jn(Jnz r)]^2 r dr)
   -------------------------------------------------------------
*/
void nrmJNser(nrmJN,n,nz,JNz)
   int n,nz;
   double *nrmJN,*JNz;
{
  int i;
  double jNp1v;

  for (i=1;i<=nz;i++) {
    jNp1v = jn(n+1,JNz[i]);
    nrmJN[i] = 2.0/sqr(jNp1v);
  }
}

/* ---------------------------------------------------------------------------
   J0trans: Feb 13, 2004
   
   Given tabulated values of J0(z_i r_j), zeros J0z, and norms nrmJ0:
     fnr given (r=0..1, j=0..nr); integrate to find Bessel coefs fnbc
 
   Note: if actual domain is r=0..R, then ensure nrmJ0 divided by R^2
   ---------------------------------------------------------------------------
*/
void j0trans(fnr,fnbc,J0array,nr,nz,J0z,nrmJ0)
   int nr,nz;
   double *fnr,*fnbc,*J0z,*nrmJ0;
   double **J0array;
{
  int i,j;
  double r,dr;

  dr = 1.0/((double) nr);

  for (i=1;i<=nz;i++) fnbc[i]=0.0;

  r = 0.5*dr;
  for (j=1;j<=nr;j++) {
    for (i=1;i<=nz;i++) 
      fnbc[i] += 0.5*(fnr[j-1]*J0array[i][j-1]+fnr[j]*J0array[i][j])*r*dr; 
    r += dr;
  }
  for (i=1;i<=nz;i++) fnbc[i] *= nrmJ0[i];
}

/* ---------------------------------------------------------------------------
   J1trans: Feb 18, 2004
   
   Given tabulated values of J1(z_i r_j), zeros J1z, and norms nrmJ1:
     fnr given (r=0..1, j=0..nr); integrate to find Bessel coefs fnbc
 
   Note: if actual domain is r=0..R, then ensure nrmJ1 divided by R^2
   ---------------------------------------------------------------------------
*/
void j1trans(fnr,fnbc,J1array,nr,nz,J1z,nrmJ1)
   int nr,nz;
   double *fnr,*fnbc,*J1z,*nrmJ1;
   double **J1array;
{
  int i,j;
  double r,dr;

  dr = 1.0/((double) nr);

  for (i=1;i<=nz;i++) fnbc[i]=0.0;

  r = 0.5*dr;
  for (j=1;j<=nr;j++) {
    for (i=1;i<=nz;i++) 
      fnbc[i] += 0.5*(fnr[j-1]*J1array[i][j-1]+fnr[j]*J1array[i][j])*r*dr; 
    r += dr;
  }
  for (i=1;i<=nz;i++) fnbc[i] *= nrmJ1[i];
}

/* ---------------------------------------------------------------------------
   jtrans: June 9, 2012
   
   Given tabulated values of JN(z_i r_j), zeros JNz, and norms nrmJN:
     fnr given (r=0..1, j=0..nr); integrate to find Bessel coefs fnbc
 
   Note: if actual domain is r=0..R, then ensure nrmJN divided by R^2
   Also note: this is identical to j0trans and j1trans, but gave name 
              jtrans for clarity
   ---------------------------------------------------------------------------
*/
void jtrans(fnr,fnbc,JNarray,nr,nz,JNz,nrmJN)
   int nr,nz;
   double *fnr,*fnbc,*JNz,*nrmJN;
   double **JNarray;
{
  int i,j;
  double r,dr;

  dr = 1.0/((double) nr);

  for (i=1;i<=nz;i++) fnbc[i]=0.0;

  r = 0.5*dr;
  for (j=1;j<=nr;j++) {
    for (i=1;i<=nz;i++) 
      fnbc[i] += 0.5*(fnr[j-1]*JNarray[i][j-1]+fnr[j]*JNarray[i][j])*r*dr; 
    r += dr;
  }
  for (i=1;i<=nz;i++) fnbc[i] *= nrmJN[i];
}

/* ---------------------------------------------------------------------------
   J0invtrans: Feb 13, 2004
   
   Given tabulated values of J0(z_i r_j):
     fnbc given (i=1..nz); integrate to find real function fnr
   ---------------------------------------------------------------------------
*/
void j0invtrans(fnr,fnbc,J0array,nr,nz)
   int nr,nz;
   double *fnr,*fnbc;
   double **J0array;
{
  int i,j;

  for (j=0;j<=nr;j++) fnr[j]=0.0;
  for (i=1;i<=nz;i++) 
  for (j=0;j<=nr;j++) fnr[j] += fnbc[i]*J0array[i][j];
}

/* ---------------------------------------------------------------------------
   J1invtrans: Feb 18, 2004
   
   Given tabulated values of J1(z_i r_j):
     fnbc given (i=1..nz); integrate to find real function fnr
   ---------------------------------------------------------------------------
*/
void j1invtrans(fnr,fnbc,J1array,nr,nz)
   int nr,nz;
   double *fnr,*fnbc;
   double **J1array;
{
  int i,j;

  for (j=0;j<=nr;j++) fnr[j]=0.0;
  for (i=1;i<=nz;i++) 
  for (j=0;j<=nr;j++) fnr[j] += fnbc[i]*J1array[i][j];
}

/* ---------------------------------------------------------------------------
   Jinvtrans: Jun 9, 2012
   
   Given tabulated values of JN(z_i r_j):
     fnbc given (i=1..nz); integrate to find real function fnr
   Note: this routine is identical to j0invtrans and j1invtrans
   ---------------------------------------------------------------------------
*/
void jinvtrans(fnr,fnbc,JNarray,nr,nz)
   int nr,nz;
   double *fnr,*fnbc;
   double **JNarray;
{
  int i,j;

  for (j=0;j<=nr;j++) fnr[j]=0.0;
  for (i=1;i<=nz;i++) 
  for (j=0;j<=nr;j++) fnr[j] += fnbc[i]*JNarray[i][j];
}

/* ---------------------------------------------------------------------------
   J0CFSERIES: Feb 13, 2004
   
   Determine J0 Bessel series coefficients from function defined for 
   r=0..R (as i=0..nr)
   
   ---------------------------------------------------------------------------
*/
void j0cfseries(fnr,nr,fnbc,nz)
  int nr,nz;
  double *fnr,*fnbc;
{
  int i;
  double *J0z,*nrmJ0;
  double **J0array;

  /* allocate space */
  J0z = dvector(1,nz);
  nrmJ0 = dvector(1,nz);
  J0array = dmatrix2(1,nz,0,nr);

  /* compute zeroes */
  for (i=1;i<=nz;i++) J0z[i]=bessj0zero(i);

  /* compute coef normalization.  Addition 1/R^2 factor then put in */
  nrmJ0ser(nrmJ0,nz,J0z);

  /* compute bessel functions with different zeroes */
  defj0array(J0array,nr,nz,J0z);
  j0trans(fnr,fnbc,J0array,nr,nz,J0z,nrmJ0);

  /* free space */
  free_dvector(J0z,1,nz);
  free_dvector(nrmJ0,1,nz);
  free_dmatrix2(J0array,1,nz,0,nr);
}

/* ---------------------------------------------------------------------------
   J1CFSERIES: Feb 18, 2004
   
   Determine J1 Bessel series coefficients from function defined for 
   r=0..R (as i=0..nr)
   
   ---------------------------------------------------------------------------
*/
void j1cfseries(fnr,nr,fnbc,nz)
  int nr,nz;
  double *fnr,*fnbc;
{
  int i;
  double *J1z,*nrmJ1;
  double **J1array;

  /* allocate space */
  J1z = dvector(1,nz);
  nrmJ1 = dvector(1,nz);
  J1array = dmatrix2(1,nz,0,nr);

  /* compute zeroes */
  for (i=1;i<=nz;i++) J1z[i]=bessj1zero(i);

  /* compute coef normalization.  Addition 1/R^2 factor then put in */
  nrmJ1ser(nrmJ1,nz,J1z);

  /* compute bessel functions with different zeroes */
  defj1array(J1array,nr,nz,J1z);
  j1trans(fnr,fnbc,J1array,nr,nz,J1z,nrmJ1);

  /* free space */
  free_dvector(J1z,1,nz);
  free_dvector(nrmJ1,1,nz);
  free_dmatrix2(J1array,1,nz,0,nr);
}

/* ---------------------------------------------------------------------------
   J0SUMSER: Feb 13, 2004
   
   Given J0 Bessel series coefficients fnbc, sum up series to get fnr
   ---------------------------------------------------------------------------
*/
void j0sumseries(fnbc,nz,fnr,nr)
  int nr,nz;
  double *fnr,*fnbc;
{
  int i;
  double *J0z;
  double **J0array;

  /* allocate space */
  J0z = dvector(1,nz);
  J0array = dmatrix2(1,nz,0,nr);

  /* define zeros and Bessel functions */
  for (i=1;i<=nz;i++) J0z[i]=bessj0zero(i);
  defj0array(J0array,nr,nz,J0z);

  /* sum up series */
  j0invtrans(fnr,fnbc,J0array,nr,nz);
 
  /* free space */
  free_dvector(J0z,1,nz);
  free_dmatrix2(J0array,1,nz,0,nr);
}

/* ---------------------------------------------------------------------------
   J1SUMSER: Feb 18, 2004
   
   Given J1 Bessel series coefficients fnbc, sum up series to get fnr
   ---------------------------------------------------------------------------
*/
void j1sumseries(fnbc,nz,fnr,nr)
  int nr,nz;
  double *fnr,*fnbc;
{
  int i;
  double *J1z;
  double **J1array;

  /* allocate space */
  J1z = dvector(1,nz);
  J1array = dmatrix2(1,nz,0,nr);

  /* define zeros and Bessel functions */
  for (i=1;i<=nz;i++) J1z[i]=bessj1zero(i);
  defj1array(J1array,nr,nz,J1z);

  /* sum up series */
  j1invtrans(fnr,fnbc,J1array,nr,nz);
 
  /* free space */
  free_dvector(J1z,1,nz);
  free_dmatrix2(J1array,1,nz,0,nr);
}

