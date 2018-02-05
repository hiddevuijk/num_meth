#ifndef GUARD_ORIENTATION_H
#define GUARD_ORIENTATION_H


#include "matrix_routines.h"
#include "lu_decomp.h"

#include <vector>
#include <math.h>

namespace orientation_f {
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
	double f3y(double y, double q, double v0,double dy)
	{
		double a = -2*B(y,q)*Bp(y,q);
		a /= 1+B(y,q)*B(y,q);
		return dy*dy*v0*a/3.;
	}
	double f3x(double y, double q, double v0,double dy)
	{
		double a = (1-B(y,q)*B(y,q))*Bp(y,q);
		a /= 1+B(y,q)*B(y,q);
		return dy*dy*v0*a/3.;
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
		
		m(N-1,0) = 1;
		m(N-1,1) = -1;
		m(N-1,N-2) = -1;
		m(N-1,N-1) = 1;	
		return m;
	}

	std::vector<double> fdy_vec(double q,double v0,double dy,int N)
	{
		std::vector<double> b(N,0);
		for(int i=1;i<N-1;++i) 
			b[i] = f3y(i*dy,q,v0,dy);
		b[0] = 0;
		b[N-1] = 0;
		return b;
	} 

	std::vector<double> fdx_vec(double q,double v0,double dy,int N)
	{
		std::vector<double> b(N,0);
		for(int i=1;i<N-1;++i) 
			b[i] = f3x(i*dy,q,v0,dy);
		b[0] = 0;
		b[N-1] = 0;
		return b;
	} 

}

double orientation(std::vector<double>& px,
	std::vector<double>& py,
	int N, double q, double v0, double Dr)
{
	double yN = 10; // period
	double dy = yN/(N-1); // interval size
	std::vector<double> b; // rhs
	std::vector<double> bb; // rhs calculated from the result
	double ex = 0; // error in px
	double ey = 0; // error in py

	Matd m = orientation_f::fd_mat(q,Dr,v0,dy,N);
	LU lu(m);

	//check determinant
	//fabs(lu.det());

	// px
	b = orientation_f::fdx_vec(q,v0,dy,N);
	px = lu.solve(b);
	bb = multiply(m,px);

	for(int i=0;i<N;++i)
		ex += (b[i]-bb[i])*(b[i]-bb[i]);

	// py
	b = orientation_f::fdy_vec(q,v0,dy,N);
	py = lu.solve(b);
	bb = multiply(m,py);
	for(int i=0;i<N;++i)
		ey += (b[i]-bb[i])*(b[i]-bb[i]);



	// return largest error
	if(ex>ey) return ex;
	else return ey;
}

#endif
