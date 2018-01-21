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
	double& operator() (std::vector<double>::size_type i,
		std::vector<double>::size_type j)
		{ return M[i][j];}
	double operator() (std::vector<double>::size_type i,
		std::vector<double>::size_type j)
		const { return M[i][j];}
	// transpose, implement by bool, and switch in i,j in operator()
	void T();
	// determinant
	double D() const;

	std::ostream& print(std::ostream&);
protected:
	std::vector<std::vector<double > > M;
	// rows, columns
	size_type r,c;

};

class Rvec: public Matd {
public:
	Rvec(int cc):Matd(1,cc) {}
	Rvec(int cc, double val): Matd(1,cc,val) {}
	std::vector<double>::size_type size() const {return c;}
	double& operator() (std::vector<double>::size_type i)
		{ return M[0][i];}
	double operator() (std::vector<double>::size_type i)
		const { return M[0][i];}
	double& operator[] (std::vector<double>::size_type i)
		{ return M[0][i];}
	double operator[] (std::vector<double>::size_type i)
		const { return M[0][i];}
};

class Cvec: public Matd {
public:
	Cvec(int rr):Matd(rr,1) {}
	Cvec(int rr, double val): Matd(rr,1,val) {}
	std::vector<double>::size_type size() const {return r;}
	double& operator() (std::vector<double>::size_type i)
		{ return M[i][0];}
	double operator() (std::vector<double>::size_type i)
		const { return M[i][0];}
	double& operator[] (std::vector<double>::size_type i)
		{ return M[i][0];}
	double operator[] (std::vector<double>::size_type i)
		const { return M[i][0];}


};

std::ostream& Matd::print(std::ostream& out)
{
	for(size_type ri=0;ri<r;++ri){
		for(size_type ci=0;ci<c;++ci){
			out << M[ri][ci];
			if(ci+1!=c) out << '\t';
		}
		out << '\n';
	}
	
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



Matd multiply(Matd a, Matd b) 
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












#endif
