/* ======================================================================
   CONVERT_XY_RT: Oct 19/95

   Convert double precision field defined on evenly spaced Cartesian grid
   to field defined on evenly spaced grid in polar co-ordinates.  Angle, t,
   is presumably given in radians.
   ======================================================================
*/
#include <stdio.h>
#include <math.h>
#include "complex.h"
#include "alloc_space.h"
#include "differentiate.h"
#include "interpolate.h"
#include "constants.h"
#include "macros.h"

void convert_xy_rt(fxy,ox,nx,oy,ny,xmin,xmax,ymin,ymax, 
                   frt,or,nr,ot,nt,rmin,rmax,tmin,tmax)
   int ox,nx,oy,ny,or,nr,ot,nt;
   double xmin,xmax,ymin,ymax,rmin,rmax,tmin,tmax;
   double **fxy,**frt;
{
   int ir,it,ix,iy,ix0,iy0;
   double xx,yy,dx,dy,dr,dt,f1,f2,f3,if1,if2,if3;
   double *x,*y,*r,*t;

   x = dvector(ox,nx);
   y = dvector(oy,ny);
   r = dvector(or,nr);
   t = dvector(ot,nt);
   
   dx = (xmax-xmin)/((double) (nx-ox));
   dy = (ymax-ymin)/((double) (ny-oy));
   dr = (rmax-rmin)/((double) (nr-or));
   dt = (tmax-tmin)/((double) (nt-ot));
   for (ix=ox;ix<=nx;ix++) x[ix] = xmin + (ix-ox)*dx;
   for (iy=oy;iy<=ny;iy++) y[iy] = ymin + (iy-oy)*dy;
   for (ir=or;ir<=nr;ir++) r[ir] = rmin + (ir-or)*dr;
   for (it=ot;it<=nt;it++) t[it] = tmin + (it-ot)*dt;

   for (ir=or;ir<=nr;ir++) 
   for (it=ot;it<=nt;it++) {
      xx = r[ir]*cos(t[it]);
      yy = r[ir]*sin(t[it]);

      ix0 = ox + (int) ((xx-xmin)/dx+0.5);
      ix0 = (ix0<ox ? ox : (ix0>nx ? nx : ix0));
      iy0 = oy + (int) ((yy-ymin)/dy+0.5);
      iy0 = (iy0<oy ? oy : (iy0>ny ? ny : iy0));

      /* determine interpolated x values at three points for different y */
      /* first four conditions extrapolate from corners of xy domain */
      if (ix0==ox && iy0==oy) {
         frt[ir][it] = fxy[ox][oy];
      } else if (ix0==ox && iy0==ny) {
         frt[ir][it] = fxy[ox][ny];
      } else if (ix0==nx && iy0==oy) {
         frt[ir][it] = fxy[nx][oy];
      } else if (ix0==nx && iy0==ny) {
         frt[ir][it] = fxy[nx][ny];
      } else if (ix0==ox) {  /* extrapolate horizontally left of xy domain */
         f1 = fxy[ox][iy0-1]; f2 = fxy[ox][iy0]; f3 = fxy[ox][iy0+1];
         frt[ir][it] = interpolate(yy, y[iy0-1],y[iy0],y[iy0+1],f1,f2,f3);
      } else if (ix0==nx) {  /* extrapolate horizontally right of xy domain */
         f1 = fxy[nx][iy0-1]; f2 = fxy[nx][iy0]; f3 = fxy[nx][iy0+1];
         frt[ir][it] = interpolate(yy, y[iy0-1],y[iy0],y[iy0+1],f1,f2,f3);
      } else if (iy0==oy) {  /* extrapolate below xy domain */
         f1 = fxy[ix0-1][oy]; f2 = fxy[ix0][oy]; f3 = fxy[ix0+1][oy];
         frt[ir][it] = interpolate(xx, x[ix0-1],x[ix0],x[ix0+1],f1,f2,f3);
      } else if (iy0==ny) {  /* extrapolate above xy domain */
         f1 = fxy[ix0-1][ny]; f2 = fxy[ix0][ny]; f3 = fxy[ix0+1][ny];
         frt[ir][it] = interpolate(xx, x[ix0-1],x[ix0],x[ix0+1],f1,f2,f3);
      } else {  /* interpolate within xy domain */
         f1 = fxy[ix0-1][iy0-1]; f2 = fxy[ix0][iy0-1]; f3 = fxy[ix0+1][iy0-1];
         if1 = interpolate(xx, x[ix0-1],x[ix0],x[ix0+1],f1,f2,f3);
         f1 = fxy[ix0-1][iy0]; f2 = fxy[ix0][iy0]; f3 = fxy[ix0+1][iy0];
         if2 = interpolate(xx, x[ix0-1],x[ix0],x[ix0+1],f1,f2,f3);
         f1 = fxy[ix0-1][iy0+1]; f2 = fxy[ix0][iy0+1]; f3 = fxy[ix0+1][iy0+1];
         if3 = interpolate(xx, x[ix0-1],x[ix0],x[ix0+1],f1,f2,f3);
         frt[ir][it] = interpolate(yy, y[iy0-1],y[iy0],y[iy0+1],if1,if2,if3);
      }
   }

   free_dvector(x,ox,nx);
   free_dvector(y,oy,ny);
   free_dvector(r,or,nr);
   free_dvector(t,ot,nt);
}

void convert_rt_xy(frt,or,nr,ot,nt,rmin,rmax,tmin,tmax,
                   fxy,ox,nx,oy,ny,xmin,xmax,ymin,ymax)
                   
   int ox,nx,oy,ny,or,nr,ot,nt;
   double xmin,xmax,ymin,ymax,rmin,rmax,tmin,tmax;
   double **fxy,**frt;
{
   int bcflg;
   int ir,it,ix,iy,ir0,it0;
   double rr,tt,dx,dy,dr,dt,f1,f2,f3,if1,if2,if3;
   double *x,*y,*r,*t;

   x = dvector(ox,nx);
   y = dvector(oy,ny);
   r = dvector(or,nr);
   t = dvector(ot,nt);
   
   bcflg = (tmin>=0.0 ? 0 : 1); /* branch cut theta=0 (0) or theta = PI (1) */

   dx = (xmax-xmin)/((double) (nx-ox));
   dy = (ymax-ymin)/((double) (ny-oy));
   dr = (rmax-rmin)/((double) (nr-or));
   dt = (tmax-tmin)/((double) (nt-ot));
   for (ix=ox;ix<=nx;ix++) x[ix] = xmin + (ix-ox)*dx;
   for (iy=oy;iy<=ny;iy++) y[iy] = ymin + (iy-oy)*dy;
   for (ir=or;ir<=nr;ir++) r[ir] = rmin + (ir-or)*dr;
   for (it=ot;it<=nt;it++) t[it] = tmin + (it-ot)*dt;

   for (ix=ox;ix<=nx;ix++) 
   for (iy=oy;iy<=ny;iy++) {
      rr = hypot(x[ix],y[iy]);
      tt = atan2(y[iy],x[ix]);       /* angle between -PI and PI */
      if (bcflg==0)
         tt = (tt<0.0 ? (tt+PI2) : tt); /* angle between 0 and 2 PI */

      ir0 = or + (int) ((rr-rmin)/dr+0.5);
      ir0 = (ir0<or ? or : (ir0>nr ? nr : ir0));
      it0 = ot + (int) ((tt-tmin)/dt+0.5);
      it0 = (it0<ot ? ot : (it0>nt ? nt : it0));

      /* determine interpolated r values at three points for different t */
      /* first four conditions extrapolate from corners of rt domain */
      if (ir0==or && it0==ot) {
         fxy[ix][iy] = frt[or][ot];
      } else if (ir0==or && it0==nt) {
         fxy[ix][iy] = frt[or][nt];
      } else if (ir0==nr && it0==ot) {
         fxy[ix][iy] = frt[nr][ot];
      } else if (ir0==nr && it0==nt) {
         fxy[ix][iy] = frt[nr][nt];
      } else if (ir0==or) {  /* extrapolate horizontally left of xy domain */
         f1 = frt[or][it0-1]; f2 = frt[or][it0]; f3 = frt[or][it0+1];
         fxy[ix][iy] = interpolate(tt, t[it0-1],t[it0],t[it0+1],f1,f2,f3);
      } else if (ir0==nr) {  /* extrapolate horizontally right of xy domain */
         f1 = frt[nr][it0-1]; f2 = frt[nr][it0]; f3 = frt[nr][it0+1];
         fxy[ix][iy] = interpolate(tt, t[it0-1],t[it0],t[it0+1],f1,f2,f3);
      } else if (it0==ot) {  /* extrapolate below xy domain */
         f1 = frt[ir0-1][ot]; f2 = frt[ir0][ot]; f3 = frt[ir0+1][ot];
         fxy[ix][iy] = interpolate(rr, r[ir0-1],r[ir0],r[ir0+1],f1,f2,f3);
      } else if (it0==nt) {  /* extrapolate above xy domain */
         f1 = frt[ir0-1][nt]; f2 = frt[ir0][nt]; f3 = frt[ir0+1][nt];
         fxy[ix][iy] = interpolate(rr, r[ir0-1],r[ir0],r[ir0+1],f1,f2,f3);
      } else {  /* interpolate within xy domain */
         f1 = frt[ir0-1][it0-1]; f2 = frt[ir0][it0-1]; f3 = frt[ir0+1][it0-1];
         if1 = interpolate(rr, r[ir0-1],r[ir0],r[ir0+1],f1,f2,f3);
         f1 = frt[ir0-1][it0]; f2 = frt[ir0][it0]; f3 = frt[ir0+1][it0];
         if2 = interpolate(rr, r[ir0-1],r[ir0],r[ir0+1],f1,f2,f3);
         f1 = frt[ir0-1][it0+1]; f2 = frt[ir0][it0+1]; f3 = frt[ir0+1][it0+1];
         if3 = interpolate(rr, r[ir0-1],r[ir0],r[ir0+1],f1,f2,f3);
         fxy[ix][iy] = interpolate(tt, t[it0-1],t[it0],t[it0+1],if1,if2,if3);
      }
   }

   free_dvector(x,ox,nx);
   free_dvector(y,oy,ny);
   free_dvector(r,or,nr);
   free_dvector(t,ot,nt);
}
