#ifndef GUARD_BISECT_H
#define GUARD_BISECT_H

#include <math.h>
#include <iostream>
#include <iomanip>
double bisect(double (*f)(double),
	double xl,double xu,
	double es, double imax)
{
	// iter = number of iterations
	int iter = 0;
	// ea = error estimate
	double ea;
	// test = f(xl)*f(xu), to check sign
	double test;
	double xr,xr_old;
	double fl = f(xl);
	double fr;
	do {

		xr_old = xr;
		xr = 0.5*(xl+xu);
		fr = f(xr);
	
		if(xr!= 0)
			ea = fabs((xr-xr_old)/xr);

		test = fl*fr;
		if(test < 0) xu = xr;
		else if(test > 0){
			xl = xr;
			fl = fr;
		} else ea = 0;
		++iter;
	} while(ea>es and iter < imax);
	return xr;
}

#endif
