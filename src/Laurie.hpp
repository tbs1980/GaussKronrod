#ifndef LAURIE_HPP
#define LAURIE_HPP

namespace gauss_konrad {

	template<typename _RealType>
	class Laurie
	{
	public:
		typedef _RealType RealType;
		typedef typename Eigen::Matrix<RealType,Eigen::Dynamic,1>  VectorType;
		typedef typename VectorType::Index IndexType;

		static void r_jacobi(const IndexType N,const RealType a,const RealType b,VectorType & alpha, VectorType & beta)
		{
			assert(a > -1);
			assert(b > -1);

			assert(alpha.rows()==beta.rows());
			assert(alpha.rows() > 0);

			assert(N>0)
			assert(alpha.rows() > floor(3*N/2))
			assert(beta.rows() > ceil(3*N/2))

			alpha(0) = (b-a)/(a+b+2.);
			beta(0) = pow(2.,(a+b+1.))*gamma(a+1.)*gamma(b+1.)/gamma(a+b+2.);

			for(IndexType n=1;n<N;++n)
			{
				RealType nab = 2.*n+a+b;
				alpha(n) = (b*b - a*a)/(nab*(nab+2.));
				beta(n) =  4.*(n+a)*(n+b)*n*(n+a+b)/(nab*nab*(nab+1.)*(nab-1.));
			}
		}

		static void r_jacobi_01(const IndexType N,const RealType a,const RealType b,VectorType & alpha, VectorType & beta)
		{
			assert(a > -1);
			assert(b > -1);

			assert(alpha.rows()==beta.rows());
			assert(alpha.rows() > 0);

			assert(N>0)
			assert(alpha.rows() > floor(3*N/2))
			assert(beta.rows() > ceil(3*N/2))

			r_jacobi(N, a, b, alpha, beta);

			for(IndexType n=0;n<N;++n)
			{
				alpha(n)  = (1+alpha(n))/2.;
			}

			beta(0) = beta(0)/pow(2,a+b+1.);

			for(IndexType n=1;n<N;++n)
			{
				beta(n) = beta(n)/4.;
			}

		}

		static void r_kronrod(const IndexType N,VectorType & a, VectorType & b)
		{
			assert(a.rows()==b.rows());
			assert(a.rows() > 0);
			assert(N>0);
			assert(2*N+1 == a.rows());

			VectorType s=VectorType::Zero(floor(N/2)+2);
			VectorType t=VectorType::Zero(floor(N/2)+2);

			t(1)=b(N+1);

			for(IndexType m=0;m<=N-2;++m)
			{
				RealType u=0;
				for(IndexType k=floor((m+1)/2);k>=0;--k)
				{
					k += 1;
					IndexType l=m-k;
					u += ( a(k+N+1) - a(l) )*t(k) + b(k+N+1)*s(k-1) - b(l)*s(k);
					s(k) = u;
				}

				VectorType swap(s);
				s=t;
				t=swap;
			}

			for(IndexType j=floor(N/2);j>=0;--j)
			{
				j += 1;
				s(j) = s(j-1);
			}

			for(IndexType m=N-1;m<=2*N-3;++m)
			{
				RealType u=0;
				IndexType j;
				for(IndexType k=m+1-N;k<=floor((m-1)/2);++k)
				{
					k += 1;
					IndexType l = m-k;
					j = N-1-l;
					u -= ( a(k+N+1) -  a(l) )*t(j) - b(k+N+1)*s(j) + b(l)*s(j+1);
					s(j) = u;
				}

				if(m % 2 == 0)
				{
					k = m/2 + 1;
					a(k+N+1) = a(k) + ( s(j) -b(k+N+1)*s(j+1) )/t(j+1);
				}
				else
				{
					k = (m+1)/2 + 1;
					b(k+N+1) =  s(j)/s(j+1);
				}

				VectorType swap(s);
				s=t;
				t=swap;
			}

			a(2*N) = a(N-1) - b(2*N)*s(1)/t(1);
		}

		
	};
}

#endif //LAURIE_HPP
