#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <GaussKonrad.hpp>
#include <unsupported/Eigen/MPRealSupport>

int test_r_kronrad(void)
{
	typedef mpfr::mpreal RealType;
	RealType::set_default_prec(256);
	//typedef double RealType;
	typedef gauss_konrad::Laurie<RealType> LauriePolicy;
	typedef LauriePolicy::IndexType IndexType;
	typedef LauriePolicy::VectorType VectorType;

	const IndexType N=10;
	const RealType a_in=0.1;
	const RealType b_in=0.2;
	VectorType a_out=VectorType::Zero(2*N);
	VectorType b_out=VectorType::Zero(2*N);

	LauriePolicy::r_jacobi_01(2*N,a_in, b_in, a_out, b_out);

	VectorType a=VectorType::Zero(2*N+1);
	VectorType b=VectorType::Zero(2*N+1);

	LauriePolicy::r_kronrod(N,a_out,b_out,a,b);

	for(IndexType i=0;i<a.rows();++i)
	{
		std::cout<<i<<"\t"<<a(i)<<"\t"<<b(i)<<std::endl;
	}

	return EXIT_SUCCESS;

}

int test_kronrad(void)
{

	typedef mpfr::mpreal RealType;
	RealType::set_default_prec(256);
	//typedef double RealType;
	typedef gauss_konrad::Laurie<RealType> LauriePolicy;
	typedef LauriePolicy::IndexType IndexType;
	typedef LauriePolicy::VectorType VectorType;

	const IndexType N=10;
	const RealType a_in=0.1;
	const RealType b_in=0.2;
	VectorType a_out=VectorType::Zero(2*N);
	VectorType b_out=VectorType::Zero(2*N);

	LauriePolicy::r_jacobi_01(2*N,a_in, b_in, a_out, b_out);

	VectorType x=VectorType::Zero(2*N+1);
	VectorType w=VectorType::Zero(2*N+1);

	LauriePolicy::kronrod(N,a_out,b_out,x,w);

	for(IndexType i=0;i<x.rows();++i)
	{
		std::cout<<i<<"\t"<<x(i)<<"\t"<<w(i)<<std::endl;
	}

	return EXIT_SUCCESS;
}

int test_mpkronrad(void)
{
	typedef mpfr::mpreal RealType;
	RealType::set_default_prec(256);
	//typedef double RealType;
	typedef gauss_konrad::Laurie<RealType> LauriePolicy;
	typedef LauriePolicy::IndexType IndexType;
	typedef LauriePolicy::VectorType VectorType;

	const IndexType N=10;

	VectorType x=VectorType::Zero(2*N+1);
	VectorType w=VectorType::Zero(2*N+1);

	LauriePolicy::mpkonrad(N,x,w);

	for(IndexType i=0;i<x.rows();++i)
	{
		std::cout<<std::setprecision(25)<<x(i)<<"\t"<<w(i)<<std::endl;
	}

	return EXIT_SUCCESS;

}

int main(void)
{
	int ret=EXIT_SUCCESS;
	ret += test_r_kronrad();
	ret += test_kronrad();
	ret += test_mpkronrad();
	return ret;
}
