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
	for(int i=0;i<N;++i) {
		m(i,i) = d;
		if(i+1<N){
			m(i+1,i) = -1;
			m(i,i+1) = -1;
		}
	}
	return m;
}

vector<double> fd_vec(double f,double T0,double TN,int N)
{
	vector<double> b(N,f);
	b[0] += T0;
	b[N-1] += TN;
	return b;
} 

int main()
{
	int N = 100;//number of intervals
	double xN = 10;
	double h = 0.01;
	double Ta = 0;
	double T0 = 200;
	double TN = 20;

	int n = N-1;
	double dx = xN/N;
	double d = 2+h*dx*dx;
	double f = h*dx*dx*Ta;
	
	Matd m = fd_mat(d,n);
	vector<double> b = fd_vec(f,T0,TN,n);
	vector<double> T(N+1),x(N+1,0);
	LU lu(m);
	vector<double> Ttemp = lu.solve(b);
	for(int i=0;i<n;++i){
		T[i+1] = Ttemp[i];
		x[i+1] = dx*(i+1);
	}
	x[N] = xN;
	T[0] = T0;
	T[N] = TN;
	ofstream Tout("T.dat");
	print_vec(T,Tout);
	ofstream xout("x.dat");
	print_vec(x,xout);
	return 0;
}


