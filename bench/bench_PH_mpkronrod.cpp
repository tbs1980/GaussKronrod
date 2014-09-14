#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <GaussKonrad.hpp>
#include <unsupported/Eigen/MPRealSupport>

// to reproduce results from
// http://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/
// see section 4. Tabulated Gauss-Kronrod weights and abscissae
// the algorithm was proposed by Laurie .
// D. P. Laurie (1997). Calculation of Gauss-Kronrod Quadrature Rules.
// Mathematics of Computation, 66(219), 1133-1145.

int main(void)
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

	std::cout<<std::fixed;
	for(IndexType i=0;i<x.rows();++i)
	{
		std::cout<<std::setprecision(25)<<x(i)<<"\t"<<w(i)<<std::endl;
	}

	return 0;
}
