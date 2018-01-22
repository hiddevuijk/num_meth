#ifndef GUARD_LU_DECOMP_H
#define GUARD_LU_DECOMP_H

class LUdecomp {
public:
	LUdecomp(const Matd&);

	std::ostream& print(std::ostream& os) { return LU.print(os);}

	template<class Vec>
	Vec solve(Vec);

private:
	Matd LU;
	Matd::size_type n;
};


LUdecomp::LUdecomp(const Matd& A)
:LU(A), n(A.rows())
{
	//check square?	
	double factor;
	for(Matd::size_type j=0;j<(n-1);++j) {
		for(Matd::size_type i=j+1;i<n;++i) {
			factor = LU(i,j)/LU(j,j);
			LU(i,j) = factor;
			for(Matd::size_type k=i+1;k<n;++k){
				LU(i,k) -= factor*LU(j,k);
			}
		}
	}
}

template<class Vec>			
Vec LUdecomp::solve(Vec b)
{
	//test size
	Vec x(n);
	Matd::size_type i,j;
	double sum;
	for(i=1;i<n;++i) {
		sum = b[i];
		for(j=1;j<i-1;++j) {
			sum -= LU(i,j)*b[j];
		}
		b[i] = sum;
	}
	x[n-1] = b[n-1]/LU(n-1,n-1);
	for(i=n-2;i>=0;--i) {
		sum = 0;
		for(j=i+1;j<n;++j) {
			sum += LU(i,j)*x[j];
		}
		x[i] = (b[i] - sum)/LU(i,i);
		if(i==0) break;
	}
	return x;
}



#endif
