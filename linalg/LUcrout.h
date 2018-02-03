#ifndef GUARD_LUCROUT_H
#define GUARD_LUCROUT_H


// error in 1 on L diagonal???


#include "matrix_routines.h"
#include <math.h>
#include <iostream>

void pivot(const Matd&,Matd::size_type,std::vector<int>&);

void LUcrout(const Matd& a, Matd& LU,
	std::vector<int>& p)
{
	Matd::size_type N = a.cols();
	Matd::size_type i,j,k;
	for(i=0;i<N;++i) p[i] = i;
	for(j=0;j<N;++j) {
		//pivot(a,j,p);
		for(i=0;i<=j;++i) {
			LU(p[i],j) = a(p[i],j);
			for(k=0;k<i;++k)
				LU(p[i],j) -= LU(p[i],k)*LU(p[k],j);
				
		}
		for(i=j+1;i<N;++i) {
			LU(p[i],j) = a(p[i],j);
			for(k=0;k<j;++k)
				LU(p[i],j) -= LU(p[i],k)*LU(p[k],j);
			LU(p[i],j) /= LU(p[j],j);
		}
		
	}

}

void pivot(const Matd& a,Matd::size_type j,std::vector<int>& p)
{
	
	double temp = fabs(a(j,j));
	Matd::size_type imax=j;
	Matd::size_type i;
	for(i=j+1;i<a.rows();++i) {
		if(fabs(a(p[i],j)) > temp) {
			temp = fabs(a(i,j));
			imax = i;
		}
	}
	i = p[j];
	p[j] = imax;
	p[imax] = i;
}

std::vector<double> LUsolve(const Matd& LU,
	const std::vector<double>& b)
{
	
	std::vector<double> x = b;
	double sum;
	Matd::size_type N = LU.cols();
	Matd::size_type i,j;
	for(i=0;i<N;++i) {
		


// faster ??
Matd LtimesU(const Matd& LU,std::vector<int> p)
{
	Matd A(LU.rows(),LU.cols(),0.);
	Matd::size_type i,j,k;
	for(i=0;i<LU.rows();++i){
		for(j=0;j<LU.cols();++j){
			for(k=0;k<=i and k<=j;++k) {
				if(i!=k) {
					A(i,j) += LU(p[i],k)*LU(p[k],j);
				} else{
					A(i,j) += LU(p[k],j);
				}
			}
		}
	}	
	return A;
}


#endif
