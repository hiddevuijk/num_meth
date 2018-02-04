#ifndef GUARD_MATRIX_ROUTINES_H
#define GUARD_MATRIX_ROUTINES_H

#include <vector>
#include <iostream>

// implement Cvec, Rvec
class Matd {
public:
	Matd(int rr,int cc)
		:r(rr),c(cc), M(rr,std::vector<double>(cc)) {}
	Matd(int rr,int cc,double val)
		:r(rr),c(cc), M(rr,std::vector<double>(cc,val)) {}

	typedef std::vector<double>::size_type size_type;
	size_type cols() const { return c;}
	size_type rows() const { return r;}

	// implement errors
	double& operator() (size_type i,size_type j)
			{ return M[i][j];}
	double operator() (size_type i,size_type j)
			const { return M[i][j];}

	void switch_rows(size_type i,size_type j);
	void switch_cols(size_type i,size_type j);

	std::ostream& print(std::ostream&);
protected:
	std::vector<std::vector<double > > M;
	// rows, columns
	size_type r,c;

};



void Matd::switch_rows(Matd::size_type i,
		Matd::size_type j)
{
	std::vector<double> temp = M[i];
	M[i] = M[j];
	M[j] = temp;
}

std::ostream& Matd::print(std::ostream& out)
{
	for(size_type ri=0;ri<r;++ri){
		for(size_type ci=0;ci<c;++ci){
			out << M[ri][ci];
			if(ci+1!=c) out << '\t';
		}
		out << '\n';
	}
	out << '\n';
	return out;	
}


/*
-----------------------------
	nonmember functions
-----------------------------
*/

//transpose, only for saure
Matd transpose(const Matd& m)
{
	Matd::size_type nr = m.rows();
	Matd::size_type nc = m.cols();
	Matd m2(nc,nr);

	for(Matd::size_type ci = 0;ci<nc;++ci) {
		for(Matd::size_type ri = 0;ri<nr;++ri) {
			m2(ci,ri) = m(ri,ci);
		}
	}

	return m2;
}

// implement vM, Mv
Matd multiply(const Matd& a,const Matd& b) 
{

	Matd c(a.rows(),b.cols());
	for(Matd::size_type i=0;i<a.rows();++i) {
		for(Matd::size_type j=0;j<b.cols();++j) {
			c(i,j) = 0.;
			for(Matd::size_type k=0;k<a.cols();++k)
				c(i,j) += a(i,k)*b(k,j);
		}
	}
	return c;
}


std::vector<double> multiply(const Matd& A,
	const std::vector<double>& b)
{
	// check dimensions
	std::vector<double> v(A.rows(),0.);
	std::vector<double>::size_type i,j;
	for(i=0;i<A.rows();++i) {
		for(j=0;j<A.cols();++j) {
			v[i] += A(i,j)*b[j];
		}
	}
	return v;
}


std::ostream& print_vec(const std::vector<double>& v,
	std::ostream& out)
{
	std::vector<double>::size_type i;
	for(i=0;i<v.size();++i) {
		out << v[i];
		if(i+1<v.size())
			out << '\t';
	}
	out << '\n';
	return out;
}






#endif
