#ifndef GUARD_GAUSS_ELIM_H
#define GUARD_GAUSS_ELIM_H

#include "matrix_routines.h"
#include <string>
#include <math.h>


template <class Vec>
void pivot(Matd& m, Vec& b, Matd::size_type i)
{
	//largest is the larges entry
	// in column i, starting from row i
	double largest = 0;
	// j is iterating through the column
	// l is the row index of the larges entry
	Matd::size_type j,l;
	for(j=i;j<m.rows();++j) {
		if(fabs(m(j,i))>largest) {
			largest = fabs(m(j,i));
			l = j;
		}
	}
	if(i!=l) {
		m.switch_rows(i,l);
		double temp = b[i];
		b[i] = b[l];
		b[l] = temp;
	}
}

template<class Vec>
Vec gauss_elim(Matd a,Vec b)
{
	Matd::size_type n = a.cols();
	std::vector<double> x(n);

	Matd::size_type k,i,j;
	double factor;
	// make all elems in col = k in row >k zero
	for(k=0;k<n-1;++k) {
		pivot(a,b,k);
		//if(fabs(a(k,k))<eps) throw
		for(i = k+1;i<n;++i) {
			factor = a(i,k)/a(k,k);
			for(j=0;j<n;++j) 
				a(i,j) -= factor*a(k,j);
			b[i] -= factor*b[k];
		}
	}
	
	x[n-1] = b[n-1]/a(n-1,n-1);
	for(i=n-2;i>=0;--i) {
		factor = b[i];
		for(j=i+1;j<n;++j) {
			factor -= a(i,j)*x[j];
		}
		x[i] = factor/a(i,i);
		if(i ==  0) break;
	}
		
	return x;
}








#endif
