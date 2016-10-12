#include <math.h>
#include "constants.h"
#include "macros.h"

typedef struct DCOMPLEX {double r,i;} dcomplex;
typedef struct ZOMPLEX {double re,im;} zomplex;

double CReal(a)
   dcomplex a;
{	
	return a.r;
}

double CImag(a)
   dcomplex a;
{	
	return a.i;
}

dcomplex Cadd(a,b)
   dcomplex a,b;
{	dcomplex c;
	c.r=a.r+b.r;
	c.i=a.i+b.i;
	return c;
}

dcomplex Csub(a,b)
   dcomplex a,b;
{	dcomplex c;
	c.r=a.r-b.r;
	c.i=a.i-b.i;
	return c;
}

dcomplex Cmul(a,b)
   dcomplex a,b;
{	dcomplex c;
	c.r=a.r*b.r-a.i*b.i;
	c.i=a.i*b.r+a.r*b.i;
	return c;
}

dcomplex Complex(re,im)
   double re,im;
{	dcomplex c;
	c.r=re;
	c.i=im;
	return c;
}

dcomplex Conjg(z)
   dcomplex z;
{	dcomplex c;
	c.r=z.r;
	c.i = -z.i;
	return c;
}

dcomplex Cdiv(a,b)
   dcomplex a,b;
{	dcomplex c;
	double r,den;
	if (dabs(b.r) >= dabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		c.r=(a.r+r*a.i)/den;
		c.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		c.r=(a.r*r+a.i)/den;
		c.i=(a.i*r-a.r)/den;
	}
	return c;
}

double Cabssqr(z)
   dcomplex z;
{	double ans;
	ans = (z.r)*(z.r) + (z.i)*(z.i);
	return ans;
}

double Cabs(z)
   dcomplex z;
{	double x,y,ans,temp;
	x=dabs(z.r);
	y=dabs(z.i);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

dcomplex Csqrt(z)
   dcomplex z;
{	dcomplex c;
	double x,y,w,r;
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r=0.0;
		c.i=0.0;
		return c;
	} else {
		x=dabs(z.r);
		y=dabs(z.i);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			c.r=w;
			c.i=z.i/(2.0*w);
		} else {
			c.i=(z.i >= 0) ? w : -w;
			c.r=z.i/(2.0*c.i);
		}
		return c;
	}
}

dcomplex Cexp(z)
   dcomplex z;
{
   dcomplex c;
   double expzr;

   expzr = exp(z.r);
   c.r = expzr*cos(z.i);
   c.i = expzr*sin(z.i);

   return c;
}

dcomplex Crpow(x, p)
   double x,p;
{
   dcomplex c;
   double xpowp;

   if (x==0.0) return( Complex(0.0,0.0) );

   xpowp = exp( p*log(dabs(x)) );

   if (x>0.0) 
      c = Complex( xpowp, 0.0 );
   else 
      c = Complex( xpowp*cos(PI*p), xpowp*sin(PI*p) );

   return c;
}

dcomplex RCmul(x,a)
   double x;
   dcomplex a;
{	dcomplex c;
	c.r=x*a.r;
	c.i=x*a.i;
	return c;
}

double Cnorm(a)
   dcomplex a;
{	return (dabs(a.r)+dabs(a.i)); }




zomplex Zadd(a,b)
   zomplex a,b;
{	zomplex c;
	c.re=a.re+b.re;
	c.im=a.im+b.im;
	return c;
}

zomplex Zsub(a,b)
   zomplex a,b;
{	zomplex c;
	c.re=a.re-b.re;
	c.im=a.im-b.im;
	return c;
}

zomplex Zmul(a,b)
   zomplex a,b;
{	zomplex c;
	c.re=a.re*b.re-a.im*b.im;
	c.im=a.im*b.re+a.re*b.im;
	return c;
}

zomplex Zomplex(re,im)
   double re,im;
{	zomplex c;
	c.re=re;
	c.im=im;
	return c;
}

zomplex Zonjg(z)
   zomplex z;
{	zomplex c;
	c.re = z.re;
	c.im = -z.im;
	return c;
}

zomplex Zdiv(a,b)
   zomplex a,b;
{	zomplex c;
	double r,den;
	if (dabs(b.re) >= dabs(b.im)) {
		r=b.im/b.re;
		den=b.re+r*b.im;
		c.re=(a.re+r*a.im)/den;
		c.im=(a.im-r*a.re)/den;
	} else {
		r=b.re/b.im;
		den=b.im+r*b.re;
		c.re=(a.re*r+a.im)/den;
		c.im=(a.im*r-a.re)/den;
	}
	return c;
}

double Zabssqr(z)
   zomplex z;
{
   return((z.re)*(z.re) + (z.im)*(z.im));
}

double Zabs(z)
   zomplex z;
{	double x,y,ans,temp;
	x=dabs(z.re);
	y=dabs(z.im);
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}

zomplex Zsqrt(z)
   zomplex z;
{	zomplex c;
	double x,y,w,r;
	if ((z.re == 0.0) && (z.im == 0.0)) {
		c.re=0.0;
		c.im=0.0;
		return c;
	} else {
		x=dabs(z.re);
		y=dabs(z.im);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.re >= 0.0) {
			c.re=w;
			c.im=z.im/(2.0*w);
		} else {
			c.im=(z.im >= 0) ? w : -w;
			c.re=z.im/(2.0*c.im);
		}
		return c;
	}
}

zomplex Zexp(z)
   zomplex z;
{
   zomplex c;
   double expzr;

   expzr = exp(z.re);
   c.re = expzr*cos(z.im);
   c.im = expzr*sin(z.im);

   return c;
}

zomplex Zrpow(x, p)
   double x,p;
{
   zomplex c;
   double xpowp;

   if (x==0.0) return( Zomplex(0.0,0.0) );

   xpowp = exp( p*log(dabs(x)) );

   if (x>0.0) 
      c = Zomplex( xpowp, 0.0 );
   else 
      c = Zomplex( xpowp*cos(PI*p), xpowp*sin(PI*p) );

   return c;
}

zomplex Zconjg(a)
   zomplex a;
{	zomplex c;
	c.re=a.re;
	c.im=-a.im;
	return c;
}

zomplex RZmul(x,a)
   double x;
   zomplex a;
{	zomplex c;
	c.re=x*a.re;
	c.im=x*a.im;
	return c;
}

double Znorm(a)
   zomplex a;
{	return (dabs(a.re)+dabs(a.im)); }
