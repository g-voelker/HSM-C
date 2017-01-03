#include <fftw3.h>
#include <math.h>
#include <stdlib.h>
#include "header.h"
#include "input.h"

void solve(fftw_plan fft, fftw_plan ifft,
           double r0, double f0, double rho0, int leap,
           double *taux, double *tauy, double *mld,
           double *freqs, fftw_complex *aux, fftw_complex *AUX){
  // declare and allocate auxiliary variables
  fftw_complex *transfer;
  double temp[2] = {0, 0};
  int nn=0;

  transfer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * DFFT_LEN);
  if (DBGFLG>2) { printf("  solveode: resetting auxiliary arrays\n"); fflush(NULL);}
  // reset forcing and initialize transfer function
  for (nn=0; nn<(365 + 62 + leap)*24; nn++){
    transfer[nn][0] = 1.0;
    transfer[nn][1] = 0.0;
    aux[nn][0] = taux[nn] / mld[nn] / rho0;
    aux[nn][1] = tauy[nn] / mld[nn] / rho0;
  }
  for (nn=(365 + 62 + leap)*24; nn<DFFT_LEN; nn++){
    transfer[nn][0] = 1.0;
    transfer[nn][1] = 0.0;
    aux[nn][0] = 0.0;
    aux[nn][1] = 0.0;
  }

  if (DBGFLG>2) { printf("  solveode: transforming in time\n"); fflush(NULL);}
  // transform in time
  fftw_execute(fft);

  if (DBGFLG>2) { printf("  solveode: setting transfer function\n"); fflush(NULL);}
  // set transfer function
  for (nn=0; nn<DFFT_LEN; nn++){
    transfer[nn][0] = r0 / (r0*r0 + (freqs[nn] + f0)*(freqs[nn] + f0));
    transfer[nn][1] = - (freqs[nn] + f0)  / (r0*r0 + (freqs[nn] + f0)*(freqs[nn] + f0));
    temp[0] = AUX[nn][0];
    temp[1] = AUX[nn][1];
    AUX[nn][0] = (temp[0] * transfer[nn][0] - temp[1] * transfer[nn][1]);
    AUX[nn][1] = (temp[1] * transfer[nn][0] + temp[0] * transfer[nn][1]);
  }

  if (DBGFLG>2) { printf("  solveode: transforming back\n"); fflush(NULL);}
  // transform back
  fftw_execute(ifft);

  for (nn=0; nn<DFFT_LEN; nn++){
    aux[nn][0] *= 1.0/DFFT_LEN;
    aux[nn][1] *= 1.0/DFFT_LEN;
  }

  fftw_free(transfer);

  if (DBGFLG>2) { printf("  solveode: return to main\n"); fflush(NULL);}
}