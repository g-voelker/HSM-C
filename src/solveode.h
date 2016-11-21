extern void solve(fftw_plan fft, fftw_plan ifft,
                  double r0, double f0,double rho0, int leap,
                  double *taux, double *tauy, double *mld,
                  double *freqs, fftw_complex *aux, fftw_complex *AUX);