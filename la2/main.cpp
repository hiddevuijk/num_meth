#include "matrix_routines.h"
#include "gauss_elim.h"
#include "lu_decomp.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

int main()
{
	Matd a(3,3);
	a(0,0) = 2;
	a(0,1) = 3;
	a(0,2) = 4;
	a(1,0) = 1;
	a(1,1) = 5;
	a(1,2) = 6;
	a(2,0) = 3;
	a(2,1) = 10;
	a(2,2) = -1;
	a.switch_rows(1,2);
	a.print(cout);
	LU lu(a);
	
	vector<double> b(3);
	b[0] = 7.85;
	b[1] = -19.3;
	b[2] = 71.4;
	vector<double> x = lu.solve(b);
	print_vec(b,cout);
	b = multiply(a,x);
	print_vec(b,cout);

	double d = lu.det();
	cout << d << endl;
	return 0;
}
