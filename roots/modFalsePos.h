#ifndef GUARD_MODFALSEPOS_H
#define GUARD_MODFALSEPOS_H

#include <math.h>

double modFalsePos(double (*f)(double),
	double xl,double xu, double es, int imax)
{
	int iter = 0;
	double xr,xr_old;
	double fl = f(xl);
	double fu = f(xu);
	double fr;
	double test;
	int iu = 0;
	int il = 0;
	double ea;
	do {
		xr_old = xr;
		xr = xu - fu*(xl-xu)/(fl-fu);
		fr = f(xr);

		if(xr != 0) ea = fabs((xr-xr_old)/xr);
		
		test = fl*fr;
		if(test<0) {
			xu = xr;
			fu = f(xu);
			iu = 0;
			++il;
			if(il>=2) fl *= 0.5;
		} else if(test > 0) {
			xl = xr;
			fl = f(xl);
			il = 0;
			++iu;
			if(iu>=2) fu *= 0.5;
		} else ea = 0.;
	} while(ea > es and iter < imax);

	return xr;
}



#endif
