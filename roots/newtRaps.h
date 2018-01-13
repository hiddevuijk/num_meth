#ifndef GUARD_NEWTRAPS_H
#define GUARD_NEWTRAPS_H

#include <math.h>
#include <iostream>

double newtRapsSave(double (*f)(double), double (*fp)(double),
	double x1,double x2, double xacc)
{
	const int imax = 100;
	double xh,xl;
	double fl = f(x1);
	double fh = f(x2);
	if(  fl*fh > 0. )
		throw("Root must be bracketed in newtRapsSave");
	if(fl == 0) return x1;
	if(fh == 0) return x2;
	if(fl < 0.0) {
		xl = x1;
		xh = x2;
	} else {
		xh = x1;
		xl = x2;
	}
	// initial guess
	double rts = .5*(x1+x2); 
	// step size before last
	double dxold = fabs(x2-x1);
	// last step size
	double dx = dxold;
	
	// f(x), f'(x)
	double fx = f(rts);
	double fpx = fp(rts);
	double temp;
	for(int i=0;i<imax;++i) {
		if((( (rts-xh)*fpx-fx)*((rts-xl)*fpx-fx) > 0.)
			|| (fabs(2.*fx) > fabs(dxold*fpx) )) {
		dxold = dx;
		dx = 0.5*(xh-xl);
		rts = xl+dx;
		if(xl == rts) return rts;
		} else {
			dxold = dx;
			dx = fx/fpx;
			temp = rts;
			rts -= dx;
			if(temp == rts) return rts;
		}
		if( fabs(dx) < xacc) return rts;
		fx = f(rts);
		fpx = fp(rts);
		if(fx<0.0) xl = rts;
		else xh = rts;
	}
	throw("max iteration exceeded in newtRapsSave");
}



double newtRapsBrac(double (*f)(double),double (*fp)(double),
	double x1, double x2, double xacc)
{
	const int imax = 20;
	double rtn = 0.5*(x1+x2);
	double fr,fpr,dx;
	for(int i=0;i<imax;++i) {
		fr = f(rtn);
		fpr = fp(rtn);
		dx = fr/fpr;
		rtn -= dx;
		if ( (x1-rtn)*(rtn-x2) < 0.)
			throw("jumped out of bracked");
		if(fabs(dx)< xacc) return rtn;
	}
	throw("Max iteration exceeded in newtRapsBrac");
}


double newtRaps(double (*f)(double),double (*fp)(double),
	double x0, double es, int imax)
{
	double xr = x0;
	double xr_old;
	double ea;
	int iter = 0;
	double fpr;
	double fpr_min = 1.e-3;
	do {
		xr_old = xr;
		fpr = fp(xr_old);
		 xr -= f(xr_old)/fpr;

		if(xr != 0) ea = fabs((xr-xr_old)/xr);
		++iter;
	} while(ea>es and iter < imax);

	return xr;
}	


double newtRapsSave2(double (*f)(double),double (*fp)(double),
	double xmin,double xmax, double xerr, int imax)
{
	double xr = 0.5*(xmin+xmax);
	xr = 1.;
	double dx;
	double dx_old;
	double xl = xmin;
	double xu = xmax;
	double fxr = f(xr);
	double fpxr = fp(xr);
	double fxl = f(xmin);
	double fxu = f(xmax);
	double e;
	for(int i=0;i<imax;++i) {
		// if f'(xr) is to small or newtRaps
		// out of range, do bisect
		dx_old = dx;
		dx_old = 1.;
		// do newtRaps
		dx = fxr/fpxr;
		xr -= dx;
		// if newtRaps steps outside bracet
		// or is not decreasing fast enough
		// do a bisection in stead
		
		if(  xr>xu || xr < xl ||
			fabs(2.*fxr) > fabs(dx_old*fpxr) ){
			dx = 0.5*(xu-xl);
			xr = xl+dx;
		}
		if(fabs(dx)<xerr) return xr; 
		fxr = f(xr);
		fpxr = fp(xr);
		if(fxr*fxl<0) {
			xu = xr;
			fxu = f(xu);
		} else {
			xl = xr;
			fxl = f(xl);
		}
		

	}
	throw("exceeded max iterations in ...");
}	


#endif
