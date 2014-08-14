#include <iostream>
#include <cstdlib>
#include <GaussKonrad.hpp>
#include <unsupported/Eigen/MPRealSupport>

int test_r_kronrad(void)
{
	//typedef mpfr::mpreal RealType;
	//RealType::set_default_prec(256);
	typedef double RealType;
	typedef gauss_konrad::Laurie<RealType> LauriePolicy;
	typedef LauriePolicy::IndexType IndexType;
	typedef LauriePolicy::VectorType VectorType;

	const IndexType N=10;
	const RealType a_in=0;
	const RealType b_in=0;
	VectorType a_out=VectorType::Zero(2*N);
	VectorType b_out=VectorType::Zero(2*N);

	LauriePolicy::r_jacobi_01(2*N,a_in, b_in, a_out, b_out);

	VectorType a=VectorType::Zero(2*N+1);
	VectorType b=VectorType::Zero(2*N+1);

	//std::cout<<floor(1/2)<<std::endl;

	LauriePolicy::r_kronrod(N,a_out,b_out,a,b);

	return EXIT_SUCCESS;

	for(IndexType i=0;i<a_out.rows();++i)
	{
		std::cout<<i<<"\t"<<a_out(i)<<"\t"<<b_out(i)<<std::endl;
	}


}

int main(void)
{
	int ret=EXIT_SUCCESS;
	ret += test_r_kronrad();
	return ret;
}