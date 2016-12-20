//
// Created by georg on 16/12/16.
//

#include "structs.h"

double davg(double *data, int ndata){
  int ii;
  double avg=0.0;
  for (ii=0; ii<ndata; ii++) avg += data[ii];
  return(avg/ndata);
}

double davg2(double **data, int n1, int n2){
  int ii, jj;
  double avg=0.0;
  for (ii=0; ii<n1; ii++) {
    for (jj=0; jj<n2; jj++) {
      avg += data[ii][jj];
    }
  }
  return(avg/(n1*n2));
}

double davg3(double ***data, int n1, int n2, int n3){
  int ii, jj, kk;
  double avg=0.0;
  for (ii=0; ii<n1; ii++) {
    for (jj=0; jj<n2; jj++) {
      for (kk=0; kk<n2; kk++) {
        avg += data[ii][jj][kk];
      }
    }
  }
  return(avg/(n1*n2*n3));
}

double dxmax(double *xx, double *yy, int index){
  int ii;
  double max=yy[0], xmax=xx[0];
  for (ii=0; ii<index; ii++){
    if (yy[ii]>max) xmax = xx[ii];
  }
  return(xmax);
}

double wavg2(dat2d *lsmask, dat3d *ww, int nt, int nxmin, int nymin){
  int ii, jj;
  double avg=0.0;
  int norm = 0;
  for (ii=1; ii<ww->nlat-1; ii++) {
    for (jj=1; jj<ww->nlon-1; jj++) {
      if (lsmask->data[ii+nymin][jj+nxmin]==0.0) {
        avg += ww->data[ii][jj][nt];
        norm++;
      }
    }
  }
  return (avg/norm);
}