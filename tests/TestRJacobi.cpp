#include <GaussKonrad.hpp>
#include <iostream>
#include <cstdlib>
#include <unsupported/Eigen/MPRealSupport>

int test_r_jacobi(void)
{
	typedef mpfr::mpreal RealType;
	RealType::set_default_prec(256);
	typedef gauss_konrad::Laurie<RealType> LauriePolicy;
	typedef LauriePolicy::IndexType IndexType;
	typedef LauriePolicy::VectorType VectorType;

	const IndexType N=10;
	const RealType a=0.5;
	const RealType b=0.7;
	VectorType alpha=VectorType::Zero(2*N+1);
	VectorType beta=VectorType::Zero(2*N+1);

	LauriePolicy::r_jacobi(3*N/2+1,a, b, alpha, beta);

	for(IndexType i=0;i<alpha.rows();++i)
	{
		std::cout<<i<<"\t"<<alpha(i)<<"\t"<<beta(i)<<std::endl;
	}

	return EXIT_SUCCESS;
}


int test_r_jacobi_01(void)
{
	typedef mpfr::mpreal RealType;
	RealType::set_default_prec(256);
	typedef gauss_konrad::Laurie<RealType> LauriePolicy;
	typedef LauriePolicy::IndexType IndexType;
	typedef LauriePolicy::VectorType VectorType;

	const IndexType N=10;
	const RealType a=0.5;
	const RealType b=0.7;
	VectorType alpha=VectorType::Zero(2*N+1);
	VectorType beta=VectorType::Zero(2*N+1);

	LauriePolicy::r_jacobi_01(3*N/2+1,a, b, alpha, beta);

	for(IndexType i=0;i<alpha.rows();++i)
	{
		std::cout<<i<<"\t"<<alpha(i)<<"\t"<<beta(i)<<std::endl;
	}

	return EXIT_SUCCESS;
}

int main(void)
{
	int ret=0;
	//ret += test_r_jacobi();
	ret += test_r_jacobi_01();
	return ret;
}
