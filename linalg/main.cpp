#include "matrix_routines.h"
#include "gauss_elim.h"

#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>


using namespace std;


int main()
{
	cout << setprecision(16);
	Matd m(2,2,0);
	m(0,0) = 0.0003;
	m(0,1) = 3;
	m(1,0) = 1;
	m(1,1) = 1;
	m.print(cout);
	cout<<endl;
	Cvec c(2);
	c(0) = 2.0001;
	c(1) = 1;
	c.print(cout);

	cout <<endl << endl;
	Cvec x = gauss_elim(m,c);
	x.print(cout) ;
	cout << endl;
	m = multiply(m,x);
	m.print(cout);
	
	return 0;
}
