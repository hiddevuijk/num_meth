#include "matrix_routines.h"
#include "orientation.h"

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <math.h>

using namespace std;



int main()
{
	int N = 5001;//number of intervals
	double q = 4.;
	double dr = 20.;
	double v0 = 10.;

	vector<double> y(N,0);
	vector<double> px,py;
	double e;
	e = orientation(px,py,N,q,v0,dr);
	cout << "error = " << e << endl;



	double dy = 10./(N-1);
	for(int i=0;i<N;++i){
		y[i] = dy*i;
	}
	ofstream pxout("px.dat");
	print_vec(px,pxout);
	ofstream pyout("py.dat");
	print_vec(py,pyout);
	ofstream yout("y.dat");
	print_vec(y,yout);
	return 0;
}


