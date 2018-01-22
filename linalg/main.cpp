#include "matrix_routines.h"
#include "gauss_elim.h"
#include "lu_decomp.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>


using namespace std;


int main()
{
	Matd A(3,3);
	A(0,0) = 3.;
	A(0,1) = -.1;
	A(0,2) = -.2;
	A(1,0) = .1;
	A(1,1) = 7;
	A(1,2) = -.3;
	A(2,0) = .3;
	A(2,1) = 0.2;
	A(2,2) = 10;


	Cvec b(3);
	b[0] = 7.85;
	b[1] = -19.3;
	b[2] = 71.4;
	cout<<setprecision(5);
	A.print(cout);	
	cout<<endl;
	LUdecomp lu(A);
	cout<<endl;

	Cvec x = lu.solve(b);
	x.print(cout);

	Matd bb = multiply(A,x);
	cout<<endl;
	bb.print(cout);

	return 0;
}
