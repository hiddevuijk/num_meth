#include "matrix_routines.h"
#include "lu_decomp.h"

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <math.h>

using namespace std;

const double pi = acos(-1);

double B(double y,double q)
{
	 return q*sin(2*pi*y/10.);
}
double Bp(double y, double q)
{
	return q*2*pi*cos(2*pi*y/10)/10;
}

double f0(double y,double q,double dy)
{

	double a = B(y,q)*Bp(y,q);
	a /= 1+B(y,q)*B(y,q);
	return 1-a*dy;
}
double f1(double y,double q,double dr,double dy)
{
	double a = 1+B(y,q)*B(y,q);
	return 2+2*dy*dy*dr*a;
}
double f2(double y,double q, double dy)
{

	double a = B(y,q)*Bp(y,q);
	a /= 1+B(y,q)*B(y,q);
	return 1+a*dy;
}
double f3y(double y, double q, double v0dr,double dy)
{
	double a = -2*B(y,q)*Bp(y,q);
	a /= 1+B(y,q)*B(y,q);
	return dy*dy*v0dr*a/6.;
}

Matd fd_mat(double q,double dr, double v0,double dy,int N)
{
	Matd m(N,N,0);
	for(int i=1;i<N-1;++i) {
		m(i,i) = -f1(i*dy,q,dr,dy);
		if(i+1<N){
			m(i+1,i) = f0(i*dy,q,dy);
			m(i,i+1) = f2(i*dy,q,dy);
		}
	}

	m(1,0) = f0(dy,q,dy);

	m(0,0) = 1;
	m(0,N-1) = -1;
	
	m(N-1,0) = -1;
	m(N-1,1) = 1;
	m(N-1,N-2) = 1;
	m(N-1,N-1) = -1;	
	return m;
}

vector<double> fd_vec(double q,double dr,double v0,double dy,int N)
{
	vector<double> b(N,0);
	for(int i=1;i<N-1;++i) 
		b[i] = f3y(i*dy,q,v0/dr,dy);
	b[0] = 0;
	b[N-1] = 0;
	return b;
} 

int main()
{
	int N = 101;//number of intervals
	double yN = 10;
	double q = 2.;
	double dr = 20.;
	double v0 = 5.;

	double dy = yN/(N-1);

	cout << setprecision(3);	
	Matd m = fd_mat(q,dr,v0,dy,N);
	//m.print(cout);
	vector<double> b = fd_vec(q,dr,v0,dy,N);
	vector<double> y(N,0);
	LU lu(m);
	//lu.print(cout);
	cout << endl;
	cout << "det = " << lu.det() << endl;
	cout << endl;
	vector<double> p = lu.solve(b);
	vector<double> bb = multiply(m,p);
	double e=0;
	for(int i=0;i<N;++i)
		e+= (b[i]-bb[i])*(b[i]-bb[i]);
	cout << "error = " << e << endl;



	for(int i=0;i<N;++i){
		y[i] = dy*i;
	}
	ofstream pout("p.dat");
	print_vec(p,pout);
	ofstream yout("y.dat");
	print_vec(y,yout);
	return 0;
}


