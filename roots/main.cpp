#include "bisect.h"
#include "modFalsePos.h"
#include "fixpt.h"
#include "newtRaps.h"

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <math.h>

using namespace std;

double f(double x)
{
	return cos(x-1.);
}

double fp(double x)
{
	return -sin(x-1.);
}

int main()
{

	double x0;

	x0 = bisect(f,-1.,0.,1.e-9,10000);
	cout << x0 << '\t' << f(x0) << endl;

	x0 = newtRapsSave2(f,fp,-1.1,1.1,1.e-9,10000);
	cout << x0 << '\t' << f(x0) << endl;


	return 0;
}


