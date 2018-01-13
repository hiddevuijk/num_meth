#ifndef GUARD_FIXPT_H
#define GUARD_FIXPT_H

#include <math.h>

double fixpt(double (*f)(double),
	double x0, double es, int imax)
{
	double xr = x0;
	double xr_old;
	double ea;
	int iter = 0;
	do {
		xr_old = xr;
		xr = f(xr_old)+xr_old;
		if(xr != 0) ea = fabs((xr-xr_old)/xr);
		++iter;
	} while(ea>es and iter < imax);

	return xr;
}	


#endif
