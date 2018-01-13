#ifndef GUARD_integrator_h
#define GUARD_integrator_h

#include <vector>
#include <math.h>

namespace intDef {
	typedef void (*DER1)(const double&, const std::vector<double>&, std::vector<double>&);
	typedef void (*SS1)(double&,std::vector<double>&,double, DER1);
	typedef void (*DER2)(const double&, const double&, double&);
	typedef void (*SS2)(double&,double&,double, DER2);

	typedef void (*adaptSS2)(const double&,const double&,const double&,
		const double&,double&,double&, DER2);


}
/*
	Driver for intagration of ODE with multiple dependent variables.

	Integration from x until xend, with stepsize h.
	Initial conditions: x,y.
	SingleStep: integration routine.
	derivs: the ODEs.

*/


void integrator(double& x, std::vector<double>& y, double h,const double& xend,
	intDef::SS1 singleStep, intDef::DER1 derivs)
{
	if(xend>x){
		while(x<xend) {
			if(xend-x<h)
				h = xend - x;
			singleStep(x,y,h,derivs);
		}
	}
	else if(xend>x) {
		while(x>xend) {
			if(xend-x>(-1*h))
				h = xend - x;
			singleStep(x,y,-1*h,derivs);
		}
	}
}


/*
	Driver for intagration of ODE with a single dependent variable.

	Integration from x until xend, with stepsize h.
	Initial conditions: x,y.
	SingleStep: integration routine.
	derivs: the ODE.

*/

void integrator(double& x,double& y, double h,const double& xend,
	intDef::SS2 singleStep, intDef::DER2 derivs)
{
	if(xend>x){
		while(x<xend) {
			if(xend-x<h)
				h = xend - x;
			singleStep(x,y,h,derivs);
		}
	}
	else if(xend>x) {
		while(x>xend) {
			if(xend-x>(-1*h))
				h = xend - x;
			singleStep(x,y,-1*h,derivs);
		}
	}
}


void adaptIntegrator(double& x, double& y, double& xend,
	intDef::adaptSS2 singleStep, intDef::DER2 derivs)
{

	static const int maxstep = 100;
	double h = 0.01;
	static const double tiny = 1.e-30;
	double eps = 5.e-8;

	double yerr;
	double emax;
	double ytemp;
	double yscal;
	int istep = 0;

	double dy;
	while( istep<maxstep and x <= xend) {
		++istep;

		derivs(x,y,dy);
		yscal=fabs(y)+fabs(h*dy) + tiny;
		if(x+h>xend) h = xend-x;

		while(true) {	
			singleStep(x,y,dy,h,ytemp,yerr,derivs);
			emax = fabs(yerr/yscal/eps);
			if(emax<=1) break;
			h = fmax(fabs(0.9*h*pow(emax,-0.25)),0.25*fabs(h));
				
		}
		x += h;
		if(emax>1.89e-4) h=0.9*h*pow(emax,-0.2);
		else h *= 4;
		y = ytemp;
	}
}	

		
		






#endif
