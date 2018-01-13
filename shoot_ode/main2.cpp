/*

main program for testing the Runge-Kutta integration scheme

*/


//#include "integrator.h"
//#include "eulerInt.h"
//#include "rkInt.h"
#include "adaptRK.h"

#include "bisect.h"

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <math.h>

using namespace::std;

void der(const double&x , const vector<double>& y, vector<double>& dydx)
{
	double h = 5.e-8;
	double Ta = 20.;
	dydx[0] = y[1];
	dydx[1] = h*(y[0]-Ta)*(y[0]-Ta)*(y[0]-Ta)*(y[0]-Ta);
}

double f(double z0)
{
	vector<double> y(2);
	double x = 0.;
	double xend = 10.;
	y[0] = 40;
	y[1] = z0;
	adaptIntegrator(x,y,xend,adaptRKInt,der);
	return y[1]-2*z0;
}


int main()
{
	double zinit = 0;
	double z0 = bisect(f,0.,0.1,1.e-14,10000);
	cout << z0<< endl << f(z0) << endl;
	return 0;
}

