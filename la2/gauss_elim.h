#ifndef GUARD_GAUSS_ELIM_H
#define GUARD_GAUSS_ELIM_H

#include "matrix_routines.h"

#include <vector>
#include <math.h>

#include <iostream>

void pivot(Matd&,std::vector<double>&,int);
void pivot(const Matd&,int,std::vector<int>&);
std::vector<double> gauss_elim(Matd a,
	std::vector<double> b)
{
	int n = b.size();
	std::vector<double> x(n,0);
	double factor;
	// get upper triangular matrix
	for(int j=0;j<n-1;++j) {
		pivot(a,b,j);
		for(int i=j+1;i<n;++i) {
			factor = a(i,j)/a(j,j);
			for(int k=j;k<n;++k) {
				a(i,k) -= factor*a(j,k);
			}
			b[i] -= factor*b[j];
		}
	}

	// back substitution
	x[n-1] = b[n-1]/a(n-1,n-1);
	for(int i=n-2;i>=0;--i) {
		x[i] = b[i];
		for(int j=n-1;j>i;--j)
			x[i] -= x[j]*a(i,j);
		x[i] /= a(i,i);
	}

	return x;
}

std::vector<double> gauss_elim2(Matd a,
	std::vector<double> b)
{
	int n = b.size();
	std::vector<double> x(n,0);
	double factor;
	std::vector<int> p(n);
	for(int i=0;i<n;++i) p[i]=i;

	// get upper triangular matrix
	for(int j=0;j<n-1;++j) {
		pivot(a,j,p);
		for(int i=j+1;i<n;++i) {
			factor = a(p[i],j)/a(p[j],j);
			for(int k=j;k<n;++k) {
				a(p[i],k) -= factor*a(p[j],k);
			}
			b[p[i]] -= factor*b[p[j]];
		}
	}

	// back substitution
	x[n-1] = b[p[n-1]]/a(p[n-1],n-1);
	for(int i=n-2;i>=0;--i) {
		x[i] = b[p[i]];
		for(int j=n-1;j>i;--j)
			x[i] -= x[j]*a(p[i],j);
		x[i] /= a(p[i],i);
	}

	return x;
}
void pivot(const Matd& a,int j,std::vector<int>& p)
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


void pivot(Matd& a,std::vector<double>& b, int j)
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
		a.switch_rows(j,imax);
		temp = b[j];
		b[j] = b[imax];
		b[imax] = temp;
	}
}










#endif
