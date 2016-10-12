#include <math.h>
#include "complex.h"
#include "alloc_space.h"
#include "linalg.h"

void usrfun(double *x, double *fn, double **jac);

#define FREERETURN {free_dmatrix2(jac,1,n,1,n);free_dvector(fn,1,n);\
	free_dvector(p,1,n);free_ivector(indx,1,n);return;}

void mnewt(int ntrial, double *x, int n, double tolx, double tolf)
{
	int k,i,*indx;
	double errx,errf,d,*bet,**alpha;

	indx=ivector(1,n);
	p=dvector(1,n);
	fn=dvector(1,n);
	jac=dmatrix2(1,n,1,n);
	for (k=1;k<=ntrial;k++) {
		usrfun(x,fn,jac);
		errf=0.0;
		for (i=1;i<=n;i++) errf += dabs(fn[i]);
		if (errf <= tolf) FREERETURN
		for (i=1;i<=n;i++) p[i] = -fn[i];
		ludcmp(jac,n,indx,&d);
		lubksb(jac,n,indx,p);
		errx=0.0;
		for (i=1;i<=n;i++) {
			errx += dabs(p[i]);
			x[i] += p[i];
		}
		if (errx <= tolx) FREERETURN
	}
	FREERETURN
}

#undef FREERETURN
