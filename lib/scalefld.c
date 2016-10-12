/* ==================================================================
   Scale real field defined on 0..n1,0..m1 grid to real field defined
   on 0..n2,0..m2 grid.  (0 and n1, 0 and m1, 0 and n2, 0 and m2 are
   identical).
   ==================================================================
 */
#include <stdio.h>
#include <math.h>
#include "complex.h" 
#include "alloc_space.h" 
#include "interpolate.h" 
#include "differentiate.h" 

#define PI 3.14159
/* =======================================================================
   Scaling is performed by averaging over nearest two points of vector.
   This method is O(z) accurate and is O(N^2) running time.
   =======================================================================
*/
void scalevecO1(f1,n1, f2,n2)
   int n1,n2;
   double *f1,*f2;
{
   int n,iin;
   double nfrac,in,dxwt;
   
   if (n1<1 || n2<1) myerror("SCALEVECO1: indices must be positive");

   if (n1==n2) {
     for (n=1;n<=n2;n++) f2[n]=f1[n];
   } else if (n1==1) {
     for (n=1;n<=n2;n++) f2[n]=f1[1];
   } else {
      if (n2==1) {
        f2[1]=f1[1]; 
      } else {
        nfrac = n1/((double) n2);
        for (n=1;n<n2;n++) {
          in = n*nfrac;
          iin = (int) in;
          dxwt = (in-iin);
          f2[n] = (1.0-dxwt)*f1[iin] + dxwt*f1[iin+1];
        }
        f2[n2]=f1[n1];
      }
   } 
}

/* =======================================================================
   Scaling is performed by cubic splines of a vector 
   =======================================================================
*/
void scalefn(x1,f1,n1, x2,f2,n2)
   int n1,n2;
   double *x1,*f1,*x2,*f2;
{
   int n;
   double yp1,ypn;
   double *f1s;
   
   if (n1<1 || n2<1) myerror("SCALEFNO1: indices must be positive");

   f1s = dvector(1,n1);

   if (n1==1) {
     for (n=1;n<=n2;n++) f2[n]=f1[1];
   } else if (n1==2) {  /* linear interpolation */
     for (n=1;n<=n2;n++)
       f2[n] = (f1[2]-f1[1])*(x2[n]-x1[1])/(x1[2]-x1[1]) + f1[1];
   } else {
     yp1 = intplderiv(x1[1], x1[1],x1[2],x1[3], f1[1],f1[2],f1[3]);
     ypn = intplderiv(x1[n1],x1[n1],x1[n1-1],x1[n1-2],f1[n1],f1[n1-1],f1[n1-2]);
     myspline(x1,f1,n1, yp1,ypn, f1s);
     for (n=1;n<=n2;n++) mysplint(x1,f1,f1s,n1, x2[n],&f2[n]);
   } 

   free_dvector(f1s,1,n1);
}

/* =======================================================================
   ScaleShortFldO1: May 03

   Scaling of matrix of short integers:
   Average over nearest four points of source grid and assign value of
   closest integer.
   This method is O(z) accurate and is O(N^2) running time.
   =======================================================================
*/
void scaleshortfldO1(f1,n1,m1, f2,n2,m2)
   int n1,m1,n2,m2;
   short **f1,**f2;
{
   int n,m,iin,iim;
   double nfrac,mfrac,in,im,dxwt,dzwt;
   double avg;
   
   if (n1<0 || n2<=0 || m1<0 || m2<=0 )
      myerror("SCALESHORTFLDO1: indices may not be negative");

   if (n1==n2 && m1==m2) {    /* matrices identical */
      for (n=0;n<=n1;n++)
      for (m=0;m<=m1;m++) f2[n][m]=f1[n][m];
   } else if (n1==0 && m1==0) {  /* input matrix has one point alone */
      for (n=0;n<=n2;n++)
      for (m=0;m<=m2;m++) f2[n][m]=f1[0][0];
   } else if (n1==0 && m1>0) {   /* input matrix has one row alone */
      if (m2==0) {
         for (n=0;n<=n2;n++) f2[n][0]=f1[0][0];
      } else if (m2==1) {
         for (n=0;n<=n2;n++) {
            f2[n][0]=f1[0][0];
            f2[n][1]=f1[0][m1];
         }
      } else {
         mfrac = ( (m1==m2) ? 1.0 : m1/(double) m2);
         for (n=0;n<=n2;n++) {
            for (m=0;m<m2;m++) {
               im = m*mfrac;
               iim = (int) im;
               dzwt = (im-iim);
               avg = (1.0-dzwt)*f1[0][iim]+dzwt*f1[0][iim+1];
               f2[n][m] = (short) avg;
            }
            f2[n][m2] = f1[n][m1];
         }
      }
   } else if (m1==0 && n1>0) {   /* input matrix has one column alone */
      if (n2==0) {
         for (m=0;m<=m2;m++) f2[0][m]=f1[0][0];
      } else if (n2==1) {
         for (m=0;m<=m2;m++) {
            f2[0][m]=f1[0][0];
            f2[1][m]=f1[n1][0];
         }
      } else {
         nfrac = ( (n1==n2) ? 1.0 : n1/(double) n2);
         for (m=0;m<=m2;m++) {
            for (n=0;n<n2;n++) {
               in = n*nfrac;
               iin = (int) in;
               dxwt = (in-iin);
               avg = (1.0-dxwt)*f1[iin][0] + dxwt*f1[iin+1][0];
               f2[n][m] = (short) avg;
            }
            f2[n2][m]=f1[n1][m];
         }
      }
   } else {
      nfrac = ( (n1==n2) ? 1.0 : n1/(double) n2 );
      mfrac = ( (m1==m2) ? 1.0 : m1/(double) m2 );
      for (n=0;n<n2;n++)
      for (m=0;m<m2;m++) {
         /* ensure in,im lie between 0..n1, 0..m1 respectively */
         in = n*nfrac;
         im = m*mfrac;
         iin = (int) in;
         iim = (int) im;
         dxwt = (in-iin);
         dzwt = (im-iim);
         avg = (1.0-dxwt)*(1.0-dzwt)*f1[iin][iim] + 
               dxwt*(1.0-dzwt)*f1[iin+1][iim] +
               (1.0-dxwt)*dzwt*f1[iin][iim+1] +
               dxwt*dzwt*f1[iin+1][iim+1];
         f2[n][m] = (short) avg;
      }
      for (m=0;m<=m2;m++) f2[n2][m]=f2[n2-1][m];
      for (n=0;n<=n2;n++) f2[n][m2]=f2[n][m2-1];
   }
}

/* =======================================================================
   Scaling is performed by averaging over nearest four points of source grid.
   This method is O(z) accurate and is O(N^2) running time.

   Subroutine assumes f1 is (0..n1, 0..m1) and f2 is (0..n2,0..m2)
   =======================================================================
*/
void scalefldO1(f1,n1,m1, f2,n2,m2)
   int n1,m1,n2,m2;
   double **f1,**f2;
{
   int n,m,iin,iim;
   double nfrac,mfrac,in,im,dxwt,dzwt;
   
   if (n1<0 || n2<0 || m1<0 || m2<0 )
      myerror("SCALEFLDO1: indices may not be negative");

   if (n1==n2 && m1==m2) {
      for (n=0;n<=n1;n++)
      for (m=0;m<=m1;m++) f2[n][m]=f1[n][m];
   } else if (n1==0 && m1==0) {
      for (n=0;n<=n2;n++)
      for (m=0;m<=m2;m++) f2[n][m]=f1[0][0];
   } else if (n1==0 && m1>0) {
      if (m2==0) 
         for (n=0;n<=n2;n++) f2[n][0]=f1[0][0];
      else if (m2==1) 
         for (n=0;n<=n2;n++) {
            f2[n][0]=f1[0][0];
            f2[n][1]=f1[0][m1];
         }
      else {
         mfrac = ( (m1==m2) ? 1.0 : m1/(double) m2);
         for (n=0;n<=n2;n++) {
            for (m=0;m<m2;m++) {
               im = m*mfrac;
               iim = (int) im;
               dzwt = (im-iim);
               f2[n][m] = (1.0-dzwt)*f1[0][iim]+dzwt*f1[0][iim+1];
            }
            f2[n][m2] = f1[n][m1];
         }
      }
   } else if (m1==0 && n1>0) {
      if (n2==0) 
         for (m=0;m<=m2;m++) f2[0][m]=f1[0][0];
      else if (n2==1) 
         for (m=0;m<=m2;m++) {
            f2[0][m]=f1[0][0];
            f2[1][m]=f1[n1][0];
         }
      else {
         nfrac = ( (n1==n2) ? 1.0 : n1/(double) n2);
         for (m=0;m<=m2;m++) {
            for (n=0;n<n2;n++) {
               in = n*nfrac;
               iin = (int) in;
               dxwt = (in-iin);
               f2[n][m] = (1.0-dxwt)*f1[iin][0] + dxwt*f1[iin+1][0];
            }
            f2[n2][m]=f1[n1][m];
         }
      }
   } else {
      nfrac = ( (n1==n2) ? 1.0 : n1/(double) n2 );
      mfrac = ( (m1==m2) ? 1.0 : m1/(double) m2 );
      for (n=0;n<n2;n++)
      for (m=0;m<m2;m++) {
         /* ensure in,im lie between 0..n1, 0..m1 respectively */
         in = n*nfrac;
         im = m*mfrac;
         iin = (int) in;
         iim = (int) im;
         dxwt = (in-iin);
         dzwt = (im-iim);
         f2[n][m] = (1.0-dxwt)*(1.0-dzwt)*f1[iin][iim] + 
                    dxwt*(1.0-dzwt)*f1[iin+1][iim] +
                    (1.0-dxwt)*dzwt*f1[iin][iim+1] +
                    dxwt*dzwt*f1[iin+1][iim+1];
      }
      for (n=0;n<n2;n++) f2[n][m2]=2*f2[n][m2-1]-f2[n][m2-2];
      for (m=0;m<=m2;m++) f2[n2][m]=2*f2[n2-1][m]-f2[n2-2][m];
   }
}

/* =======================================================================
   Scaling using cubic splines to interpolate points.  This method is O(z^3)
   accurate but is O(N^4) running time.  Note: myspline/mysplint expect
   vectors ranging from 1..N.
   =======================================================================
*/
void scalefld(f1,n1,m1,i0flg, f2,n2,m2,o0flg)
   int n1,m1,n2,m2,i0flg,o0flg;
   double **f1,**f2;
{
   int n,m,n1low,m1low,n2low,m2low;
   double slope;
   double dx1,dx2,dy1,dy2,ihfdy,ihfdx,yp1,ypN,xp1,xpN;
   double *fx2,*fy2,*x,*y;
   double **ftmp;

   if (n1<0 || n2<0 || m1<0 || m2<0 )
      myerror("SCALEFLD: indices may not be negative");

   n1low = ((i0flg==0 || i0flg==1) ? 0 : 1);
   m1low = ((i0flg==0 || i0flg==2) ? 0 : 1);
   n2low = ((o0flg==0 || o0flg==1) ? 0 : 1);
   m2low = ((o0flg==0 || o0flg==2) ? 0 : 1);

   if (n1<n1low || n2<n2low || m1<m1low || m2<m2low )
      myerror("SCALEFLD: maximum index lower than lower index");

   if (n1==n2 && m1==m2 && n1low==n2low && m1low==m2low) {
      for (n=n1low;n<=n1;n++)
      for (m=m1low;m<=m1;m++) f2[n][m]=f1[n][m];
   } else if (n1==0 && m1==0) {
      for (n=0;n<=n2;n++)
      for (m=0;m<=m2;m++) f2[n][m]=f1[0][0];
   } else {        /* scaling required */
      /* allocate space */
      x = dvector(0,n1);
      y = dvector(0,m1);
      fx2 = dvector(0,n1);
      fy2 = dvector(0,m1);
      ftmp = dmatrix2(0,m2,0,n1);

      if (m1==0) {   /* and n1>0 */
         for (n=n1low;n<=n1;n++) 
         for (m=m2low;m<=m2;m++) ftmp[m][n]=f1[n][0];
      } else if (m1==1 && m1low==0) {  /* and n1>=0 */ 
         dy2 = ( m2==m2low ? 0.0 : 1.0/(double) (m2-m2low));
         for (n=n1low;n<=n1;n++) {
            slope = f1[n][1]-f1[n][0]; 
            for (m=m2low;m<=m2;m++) ftmp[m][n]=f1[n][0] + slope*m*dy2;
         }
      } else if (m1==1 && m1low==1) {  /* and n1>=0 */ 
         for (n=n1low;n<=n1;n++) 
         for (m=m2low;m<=m2;m++) ftmp[m][n]=f1[n][m1low];
      } else {    /* m1>=2  and n1>=0 */
         dy1 = 1.0/((double) (m1-m1low));
         dy2 = ( m2==m2low ? 0.0 : 1.0/((double) (m2-m2low)));
         ihfdy = 0.5/dy1;
         for (m=m1low;m<=m1;m++) y[m]=(m-m1low)*dy1;
 
         for (n=n1low;n<=n1;n++) {
            yp1 = ihfdy*(-3*f1[n][m1low]+4*f1[n][m1low+1]-f1[n][m1low+2]);
            ypN = ihfdy*(3*f1[n][m1]-4*f1[n][m1-1]+f1[n][m1-2]);
            myspline(y+m1low-1,f1[n]+m1low-1,m1-m1low+1,yp1,ypN, fy2+m1low-1);
            for (m=m2low;m<=m2;m++)
               mysplint(y+m1low-1,f1[n]+m1low-1,fy2+m1low-1,m1-m1low+1, 
                        (m-m2low)*dy2, &ftmp[m][n]);
         }
      }

      /* interpolate points in columns */
      if (n1==0) {  /* and m1>0 */
         for (m=m2low;m<=m2;m++)
         for (n=n2low;n<=n2;n++) f2[n][m]=ftmp[m][0];
      } else if (n1==1 && n2low==0) {   /* and m1>=0 */
         dx2 = ( n2==n2low ? 0.0 : 1.0/((double) (n2-n2low)));
         for (m=m2low;m<=m2;m++) {
            slope = ftmp[m][n1low+1]-ftmp[m][n1low];
            for (n=n2low;n<=n2;n++) f2[n][m]=ftmp[m][n1low]+slope*n*dx2;
         }
      } else {  /* n1>=2 and m1>=0 */
         dx1 = 1/((double) (n1-n1low));
         dx2 = ( n2==n2low ? 0.0 : 1.0/((double) (n2-n2low)) );
         ihfdx = 0.5/dx1;
         for (n=n1low;n<=n1;n++) x[n]=(n-n1low)*dx1;
 
         for (m=m2low;m<=m2;m++) {
            xp1 = ihfdx*(-3*ftmp[m][n1low]+4*ftmp[m][n1low+1]-ftmp[m][n1low+2]);
            xpN = ihfdx*(3*ftmp[m][n1]-4*ftmp[m][n1-1]+ftmp[m][n1-2]);

            myspline(x+n1low-1,ftmp[m]+n1low-1,n1-n1low+1,xp1,xpN, fx2+n1low-1);

            for (n=n2low;n<=n2;n++)
               mysplint(x+n1low-1,ftmp[m]+n1low-1,fx2+n1low-1,n1-n1low+1, 
                        (n-n2low)*dx2, &f2[n][m]);
          }
      }
      free_dvector(x,0,n1);
      free_dvector(y,0,m1);
      free_dvector(fx2,0,n1);
      free_dvector(fy2,0,m1);
      free_dmatrix2(ftmp,0,m2,0,n1);
   }
}

/* =======================================================================
   Scaling using cubic splines to interpolate points.  This method is O(z^3)
   accurate but is O(N^4) running time.  Note: myspline/mysplint expect
   vectors ranging from 1..N.
   =======================================================================
*/
void scalefield(f1,ox1,nx1,oy1,ny1,xmn1,xmx1,ymn1,ymx1, 
           f2,ox2,nx2,oy2,ny2,xmn2,xmx2,ymn2,ymx2,
           exflg)
   int ox1,nx1,oy1,ny1, ox2,nx2,oy2,ny2;
   int exflg;
   double xmn1,xmx1,ymn1,ymx1, xmn2,xmx2,ymn2,ymx2; 
   double **f1,**f2;
{
   int n,m;
   double dx1,dx2,dy1,dy2,yp1,ypN,xp1,xpN,xv,yv;
   double *fxs,*fys,*x,*y;
   double **ftmp;

   if (nx1<ox1 || ny1<oy1 || nx2<ox2 || ny2<oy2 )
     myerror("SCALEFIELD: indices incorrect");

   if (nx1==ox1 && ny1==oy1) {
     for (n=ox2;n<=nx2;n++)
     for (m=oy2;m<=ny2;m++) f2[n][m]=f1[ox1][oy1];
   } else if (nx1==ox1+1 || ny1==oy1+1) {
     myerror("SCALEFIELD: linear scaling not yet implemented");
   } else {        /* scaling required */
     /* allocate space */
     x = dvector(ox1,nx1);
     y = dvector(oy1,ny1);
     fxs = dvector(ox1,nx1);
     fys = dvector(oy1,ny1);
     ftmp = dmatrix2(oy2,ny2,ox1,nx1);

     dy1 = (ymx1-ymn1)/((double) (ny1-oy1));
     dy2 = (ymx2-ymn2)/((double) (ny2-oy2));
     for (m=oy1;m<=ny1;m++) y[m]=(m-oy1)*dy1 + ymn1;
 
     for (n=ox1;n<=nx1;n++) {
       yp1 = intplderiv(y[oy1], y[oy1],y[oy1+1],y[oy1+2],
                              f1[n][oy1],f1[n][oy1+1],f1[n][oy1+2]);
       ypN = intplderiv(y[ny1], y[ny1],y[ny1-1],y[ny1-2],
                              f1[n][ny1],f1[n][ny1-1],f1[n][ny1-2]);
       myspline(y+oy1-1,f1[n]+oy1-1,ny1-oy1+1,yp1,ypN, fys+oy1-1);
       for (m=oy2;m<=ny2;m++) {
          yv = ymn2+(m-oy2)*dy2;
          if (yv<ymn1 || yv>ymx1) { 
            if (exflg==0) { 
              ftmp[m][n]=0.0;
            } else if (exflg==1) {
              if (yv<ymn1) ftmp[m][n]=f1[n][oy1];
              else ftmp[m][n]=f1[n][ny1];
            } else {
              mysplint(y+oy1-1,f1[n]+oy1-1,fys+oy1-1,ny1-oy1+1,yv,&ftmp[m][n]);
            }
          } else {
            mysplint(y+oy1-1,f1[n]+oy1-1,fys+oy1-1,ny1-oy1+1, yv, &ftmp[m][n]);
          }
       }
     }


     /* interpolate points in columns */
     dx1 = (xmx1-xmn1)/((double) (nx1-ox1));
     dx2 = (xmx2-xmn2)/((double) (nx2-ox2));
     for (n=ox1;n<=nx1;n++) x[n]=(n-ox1)*dx1 + xmn1;

     for (m=oy2;m<=ny2;m++) {
       xp1 = intplderiv(x[ox1], x[ox1],x[ox1+1],x[ox1+2],
                              ftmp[m][ox1],ftmp[m][ox1+1],ftmp[m][ox1+2]);
       xpN = intplderiv(x[nx1], x[nx1],x[nx1-1],x[nx1-2],
                              ftmp[m][nx1],ftmp[m][nx1-1],ftmp[m][nx1-2]);
       myspline(x+ox1-1,ftmp[m]+ox1-1,nx1-ox1+1,xp1,xpN, fxs+ox1-1);
       for (n=ox2;n<=nx2;n++) {
         xv = (n-ox2)*dx2+xmn2;
         if (xv<xmn1 || xv>xmx1) { 
           if (exflg==0) { 
             f2[n][m]=0.0;
           } else if (exflg==1) {
             if (xv<xmn1) f2[n][m]=ftmp[m][ox1];
             else f2[n][m]=ftmp[m][nx1];
           } else {
             mysplint(x+ox1-1,ftmp[m]+ox1-1,fxs+ox1-1,nx1-ox1+1,xv,&f2[n][m]);
           }
         } else {
           mysplint(x+ox1-1,ftmp[m]+ox1-1,fxs+ox1-1,nx1-ox1+1,xv,&f2[n][m]);
         }
       }
     }

     free_dvector(x,ox1,nx1);
     free_dvector(y,oy1,ny1);
     free_dvector(fxs,ox1,nx1);
     free_dvector(fys,oy1,ny1);
     free_dmatrix2(ftmp,oy2,ny2,ox1,nx1);
   }
}

/* =======================================================================
   Rotate about axis using local averaging of 4 points to interpolate.
   =======================================================================
*/
void rotatefldO1(f1,ox,nx,oy,ny,xmn,xmx,ymn,ymx, 
           Ox,Oy,rot, f2, 
           exflg)
   int ox,nx,oy,ny;
   int exflg;
   double xmn,xmx,ymn,ymx, Ox,Oy,rot;
   double **f1,**f2;
{
   int n,m;
   int ixp,iyp;
   double Lx,Ly,dx,dy,xp,yp,fxp,fyp,dxwt,dywt;
   double *x,*y;

   /* allocate space */
   x = dvector(ox,nx);
   y = dvector(oy,ny);

   Lx = (xmx-xmn);
   Ly = (ymx-ymn);
   dx = Lx/((double) (nx-ox));
   dy = Ly/((double) (ny-oy));
   for (n=ox;n<=nx;n++) x[n]=(n-ox)*dx + xmn;
   for (m=oy;m<=ny;m++) y[m]=(m-oy)*dy + ymn;

   if (nx<ox || ny<oy ) myerror("ROTATEFIELD: indices incorrect");

   if (nx==ox && ny==oy) 
     myerror("ROTATEFIELD: field is vector ... cannot rotate");

   /* else sit and rotate! */
   /* assume new co-ords have same range and res. as orig. fld. */
   for (n=ox;n<=nx;n++) 
   for (m=oy;m<=ny;m++) {
     xp = Ox + (x[n]-Ox)*cos(rot)+(y[m]-Oy)*sin(rot); 
     yp = Oy - (x[n]-Ox)*sin(rot)+(y[m]-Oy)*cos(rot); 

     fxp = ((double) ox) + nx*(xp-xmn)/Lx;
     fyp = ((double) oy) + ny*(yp-ymn)/Ly;
     ixp = (int) fxp;
     iyp = (int) fyp;
     dxwt = fxp-ixp;
     dywt = fyp-iyp;

     if (ixp<ox) {
       if (exflg==0) {
         f2[n][m] = 0.0;
       } else { /* exflg=2 same as exflg=1 in this routine for now */
         if (iyp<oy) { f2[n][m] = f1[ox][oy]; }
         else if (iyp>=ny) { f2[n][m] = f1[ox][ny]; }
         else { f2[n][m] = (1.0-dxwt)*f1[ox][iyp] + dywt*f1[ox][iyp+1]; }
       }
     } else if (ixp>=nx) {
       if (exflg==0) {
         f2[n][m] = 0.0;
       } else { /* exflg=2 same as exflg=1 in this routine for now */
         if (iyp<oy) { f2[n][m] = f1[nx][oy]; }
         else if (iyp>=ny) { f2[n][m] = f1[nx][ny]; }
         else { f2[n][m] = (1.0-dywt)*f1[nx][iyp] + dywt*f1[nx][iyp+1]; }
       }
     } else {
       if (iyp<oy) { 
         if (exflg==0) {
           f2[n][m] = 0.0;
         } else { /* exflg=2 same as exflg=1 in this routine for now */
           f2[n][m] = (1.0-dxwt)*f1[ixp][oy] + dxwt*f1[ixp+1][oy];
         }
       } else if (iyp>=ny) { 
         if (exflg==0) {
           f2[n][m] = 0.0;
         } else { /* exflg=2 same as exflg=1 in this routine for now */
           f2[n][m] = (1.0-dxwt)*f1[ixp][ny] + dxwt*f1[ixp+1][ny];
         }
       } else {
         f2[n][m] = (1.0-dxwt)*(1.0-dywt)*f1[ixp][iyp] + 
                    dxwt*(1.0-dywt)*f1[ixp+1][iyp] +
                    (1.0-dxwt)*dywt*f1[ixp][iyp+1] +
                    dxwt*dywt*f1[ixp+1][iyp+1];
       }
     }
   }

   free_dvector(x,ox,nx);
   free_dvector(y,oy,ny);
}

/* =======================================================================
   Reflect about a line using local averaging of 4 points to interpolate.
   =======================================================================
*/
void reflectfldO1(f1,ox,nx,oy,ny,xmn,xmx,ymn,ymx, 
             x1,y1,x2,y2, f2, exflg)
   int ox,nx,oy,ny;
   int exflg;
   double xmn,xmx,ymn,ymx, x1,y1,x2,y2;
   double **f1,**f2;
{
   int n,m;
   int ixp,iyp;
   double Lx,Ly,dx,dy,xp,yp,fxp,fyp,dxwt,dywt;
   double xi,yi,ms;
   double *x,*y;

   /* allocate space */
   x = dvector(ox,nx);
   y = dvector(oy,ny);

   Lx = (xmx-xmn);
   Ly = (ymx-ymn);
   dx = Lx/((double) (nx-ox));
   dy = Ly/((double) (ny-oy));
   for (n=ox;n<=nx;n++) x[n]=(n-ox)*dx + xmn;
   for (m=oy;m<=ny;m++) y[m]=(m-oy)*dy + ymn;

   if (nx<ox || ny<oy ) myerror("REFLECTFIELD: indices incorrect");

   if (nx==ox && ny==oy) 
     myerror("REFLECTFIELD: field is vector ... cannot reflect");

   if (x1==x2 && y1==y2) {
     printf("(x1,y1)=(%10.3e,%10.3e)  (x2,y2)=(%10.3e,%10.3e)\n",x1,y1,x2,y2);
     myerror("REFLECTFIELD: cannot reflect about a point");
   }

   /* get reflective! */
   /* assume new co-ords have same range and res. as orig. fld. */
   for (n=ox;n<=nx;n++) 
   for (m=oy;m<=ny;m++) {
     if (x1==x2) { /* reflect about vertical axis */
       xp = 2*x1-x[n];
       yp = y[m];
     } else if (y1==y2) { /* reflect about horizontal axis */
       xp = x[n];
       yp = 2*y1-y[m];
     } else { /* reflect about horizontal axis */
       ms=(y2-y1)/(x2-x1);
       xi=(x[n]+x1*ms*ms + ms*(y[m]-y1))/(1.0+ms*ms);
       yi=(-(xi-x[n])+ms*y[m])/ms;
       xp = 2*xi-x[n];
       yp = 2*yi-y[m];
     }

     fxp = ((double) ox) + nx*(xp-xmn)/Lx;
     fyp = ((double) oy) + ny*(yp-ymn)/Ly;
     ixp = (int) fxp;
     iyp = (int) fyp;
     dxwt = fxp-ixp;
     dywt = fyp-iyp;

     if (ixp<ox) {
       if (exflg==0) {
         f2[n][m] = 0.0;
       } else { /* exflg=2 same as exflg=1 in this routine for now */
         if (iyp<oy) { f2[n][m] = f1[ox][oy]; }
         else if (iyp>=ny) { f2[n][m] = f1[ox][ny]; }
         else { f2[n][m] = (1.0-dxwt)*f1[ox][iyp] + dywt*f1[ox][iyp+1]; }
       }
     } else if (ixp>=nx) {
       if (exflg==0) {
         f2[n][m] = 0.0;
       } else { /* exflg=2 same as exflg=1 in this routine for now */
         if (iyp<oy) { f2[n][m] = f1[nx][oy]; }
         else if (iyp>=ny) { f2[n][m] = f1[nx][ny]; }
         else { f2[n][m] = (1.0-dywt)*f1[nx][iyp] + dywt*f1[nx][iyp+1]; }
       }
     } else {
       if (iyp<oy) { 
         if (exflg==0) {
           f2[n][m] = 0.0;
         } else { /* exflg=2 same as exflg=1 in this routine for now */
           f2[n][m] = (1.0-dxwt)*f1[ixp][oy] + dxwt*f1[ixp+1][oy];
         }
       } else if (iyp>=ny) { 
         if (exflg==0) {
           f2[n][m] = 0.0;
         } else { /* exflg=2 same as exflg=1 in this routine for now */
           f2[n][m] = (1.0-dxwt)*f1[ixp][ny] + dxwt*f1[ixp+1][ny];
         }
       } else {
         f2[n][m] = (1.0-dxwt)*(1.0-dywt)*f1[ixp][iyp] + 
                    dxwt*(1.0-dywt)*f1[ixp+1][iyp] +
                    (1.0-dxwt)*dywt*f1[ixp][iyp+1] +
                    dxwt*dywt*f1[ixp+1][iyp+1];
       }
     }
   }

   free_dvector(x,ox,nx);
   free_dvector(y,oy,ny);
}

/* =======================================================================
   Copy one field to another 
   =======================================================================
*/
void copyfield(f1,ox,nx,oy,ny, f2)
   int ox,nx,oy,ny;
   double **f1,**f2;
{
   int n,m;

   for (n=ox;n<=nx;n++) 
   for (m=oy;m<=ny;m++) f2[n][m] = f1[n][m];
}
