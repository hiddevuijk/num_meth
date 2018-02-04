#include "matrix_routines.h"
#include "lu_decomp.h"

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;

Matd fd_mat(double d,int N)
{
	Matd m(N,N,0);
	for(int i=1;i<N-1;++i) {
		m(i,i) = d;
		if(i+1<N){
			m(i+1,i) = -1;
			m(i,i+1) = -1;
		}
	}
	m(0,0)=1;
	m(1,0)=-1;
	m(0,N-1) = -1;
	
	m(N-1,0) = 1;
	m(N-1,1) = 1;
	m(N-1,N-2) = -1;
	m(N-1,N-1) = -1;	
	return m;
}

vector<double> fd_vec(double f,double T0,double TN,int N)
{
	vector<double> b(N,f);
	b[0] = T0;
	b[N-1] = TN;
	return b;
} 

int main()
{
	int N = 1001;//number of intervals
	double xN = 10;
	double h = 0.01;
	double Ta = 20;
	double T0 = 0;
	double TN = 0;

	double dx = xN/(N-1);
	double d = 2+h*dx*dx;
	double f = h*dx*dx*Ta;
	cout << setprecision(3);	
	Matd m = fd_mat(d,N);
	m.print(cout);
	vector<double> b = fd_vec(f,T0,TN,N);
	vector<double> x(N,0);
	LU lu(m);
	lu.print(cout);
	cout << endl;
	cout << lu.det() << endl;
	cout << endl;
	vector<double> T = lu.solve(b);
	vector<double> bb = multiply(m,T);
	double e=0;
	for(int i=0;i<N;++i)
		e+= (b[i]-bb[i])*(b[i]-bb[i]);
	cout << e << endl;



	for(int i=0;i<N;++i){
		x[i] = dx*i;
	}
	ofstream Tout("T.dat");
	print_vec(T,Tout);
	ofstream xout("x.dat");
	print_vec(x,xout);
	return 0;
}


