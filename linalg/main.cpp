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
	Matd m(3,3,0);
	m(0,0) = 5;
	m(0,1) = -.1;
	m(0,2) = -.2;
	m(1,0) = .1;
	m(1,1) = 7;
	m(1,2) = -.3;
	m(2,0) = .3;
	m(2,1) = -.2;
	m(2,2) = 10;
	m.print(cout);
	cout << endl;
	m.T();
	m.print(cout);
	Matd m2 = T(T(m));
	cout<< endl;
	m2.print(cout);
	return 0;
}
