//
// Created by georg on 1/2/2017.
//

#include <fftw3.h>
#include <stdlib.h>
#include <math.h>
#include "../lib/structs.h"
#include "../lib/dalloc.h"
#include "../lib/macros.h"
#include "header.h"
#include "input.h"

void hybrid(dat1d *lh, dat1d *ww, dat1d *NN, dat1d *Eout, double *freqs, double f0,
            fftw_complex *AUX, fftw_complex *aux,
            fftw_plan fft, fftw_plan ifft, int leap){
  double *window, avgN;
  int nt;

  if (DBGFLG>2) {printf("  hybrid: set up local window function\n"); fflush(NULL);}

  // get average stratification for window
  avgN = 0.0;
  for (nt=0; nt<(365 + leap)*24; nt++) avgN += NN->data[nt] / (365 + leap)*24;

  // set up window
  window = dalloc(window, (size_t) (365 + leap) * 24);
  for (nt=0; nt<(365 + leap) * 24; nt++){
    if ((freqs[nt]<=f0)|(freqs[nt]>=NN->data[nt])) {
      window[nt] = 0.0;
    } else {
      window[nt] = sqrt(sqrt(sqr(avgN) - sqr(freqs[nt])) *
                        sqrt(sqr(freqs[nt]) - sqr(f0)) /
                        dabs(freqs[nt]) / (avgN-f0));
    }
  }

  if (DBGFLG>2) {printf("  hybrid: filter vertical velocity\n"); fflush(NULL);}

  // reset auxiliary arrays
  for (nt=0; nt<(365 + leap) * 24; nt++) {
    aux[nt][0] = ww->data[nt];
    aux[nt][0] = AUX[nt][0] = AUX[nt][1] = 0.0;
  }



  // transform in time
  fftw_execute(fft);

  // multiply by filter
  for (nt=0; nt<(365 + leap) * 24; nt++) {
    AUX[nt][0] *= window[nt];
    AUX[nt][1] *= window[nt];
  }

  // transform back
  fftw_execute(ifft);

  for (nt=0; nt<(365 + leap) * 24; nt++) {
    Eout->data[nt] = RHO * (avgN - f0) * lh->data[nt] / PI2 * (aux[nt][0] * aux[nt][0] + aux[nt][1] * aux[nt][1]);
  }

  if (DBGFLG>2) {printf("  hybrid: return to main\n"); fflush(NULL);}
}