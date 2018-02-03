#ifndef GUARD_LU_DECOMP_H
#define GUARD_LU_DECOMP_H

/*
	to do
		refine
		inverse
		order of loops
*/

#include "matrix_routines.h"

#include <fstream>

#include <iostream>

void pivot_lu(const Matd&, int j,std::vector<int>&);


class LU{
public:
	LU(const Matd&);

	std::vector<double> solve(std::vector<double>) const;
	double det() const;
private:

	Matd m;
	Matd::size_type n;
	std::vector<int> p;

};

double LU::det() const
{
	double d=1;
	for(int i=0;i<n;++i)
		d*= m(p[i],i);
	return d;
}

std::vector<double> LU::solve(std::vector<double> b) const
{
	// check size
	int n = b.size();
	std::vector<double> x(n,0);
	// replace b with d
	for(int i=1;i<n;++i) 
		for(int j=0;j<i;++j)
			b[p[i]] -= m(p[i],j)*b[p[j]];
	
	x[n-1] = b[p[n-1]]/m(p[n-1],n-1);
	for(int i=n-2;i>=0;--i) {
		x[i] = b[p[i]];
		for(int j=i+1;j<n;++j) 
			x[i] -= m(p[i],j)*x[j];
		x[i]/=m(p[i],i);
	}
	return x;
}
LU::LU(const Matd& a)
	: m(a),n(a.rows()),p(a.rows())
{

	for(int i=0;i<n;++i) p[i]=i;

	double factor;

	// get upper triangular matrix
	for(int j=0;j<n-1;++j) {
		pivot(m,j,p);
		for(int i=j+1;i<n;++i) {
			factor = m(p[i],j)/m(p[j],j);
			m(p[i],j) = factor;
			for(int k=j+1;k<n;++k) {
				m(p[i],k) -= factor*m(p[j],k);
			}
		}
	}

}

void pivot_lu(const Matd& a, int j,std::vector<int>& p)
{
	double temp = fabs(a(j,j));
	int imax=j;
	for(int i=j;i<a.rows();++i) {
		if(fabs(a(i,j))>temp) {
			imax = i;
			temp = fabs(a(i,j));
		}
	}
	if(imax!=j) {
		int ti = p[imax];
		p[imax] = p[j];
		p[j] = ti;
	}
}



#endif
