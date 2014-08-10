#include <iostream>
#include <unsupported/Eigen/MPRealSupport>

int main(void)
{
	mpfr::mpreal::set_default_prec(256);

	mpfr::mpreal x=1.5;
	mpfr::mpreal gx=mpfr::gamma(x);

	std::cout<<x<<"\t"<<gx<<std::endl;

	return 0;
}