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
	Matd m(2,3,0);
	m(0,0) = 5;
	m(0,1) = -.1;
	m(0,2) = -.2;
	m(1,0) = .1;
	m(1,1) = 7;
	m(1,2) = -.3;
	cout << endl;
	m.print(cout);
	cout << endl << m.cols() << endl << m.rows() << endl <<endl;
	Matd m2 = transpose(m);
	cout<< endl;
	m2.print(cout);
	cout << endl << m2.cols() << endl << m2.rows() << endl <<endl;
	return 0;
}
