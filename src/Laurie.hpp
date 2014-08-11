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

		static void r_jacobi(const RealType a,const RealType b,VectorType & alpha, VectorType & beta)
		{
			assert(a > -1);
			assert(b > -1);

			assert(alpha.rows()==beta.rows());
			assert(alpha.rows() > 0);

			const IndexType N=alpha.rows();
			
			alpha(0) = (b-a)/(a+b+2.);
			beta(0) = pow(2.,(a+b+1.))*gamma(a+1.)*gamma(b+1.)/gamma(a+b+2.);

			for(IndexType n=1;n<N;++n)
			{
				RealType nab = 2.*n+a+b;
				alpha(n) = (b*b - a*a)/(nab*(nab+2.));
				beta(n) =  4.*(n+a)*(n+b)*n*(n+a+b)/(nab*nab*(nab+1.)*(nab-1.));
			}
		}

		static void r_jacobi_01(const RealType a,const RealType b,VectorType & alpha, VectorType & beta)
		{
			assert(a > -1);
			assert(b > -1);

			assert(alpha.rows()==beta.rows());
			assert(alpha.rows() > 0);

			const IndexType N=alpha.rows();

			r_jacobi( a, b, alpha, beta);

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
/*
function ab=r_kronrod(N,ab0)
if length(ab0)<ceil(3*N/2)+1, error('array ab0 too short'), end
a=mp(zeros(2*N+1,1)); b=a;						% P.H.
k=0:floor(3*N/2); a(k+1)=ab0(k+1,1);
k=0:ceil(3*N/2); b(k+1)=ab0(k+1,2);
s=mp(zeros(floor(N/2)+2,1)); t=s; t(2)=b(N+2);	% P.H.
for m=0:N-2,
  k=floor((m+1)/2):-1:0; l=m-k;
  s(k+2)=cumsum((a(k+N+2)-a(l+1)).*t(k+2)+b(k+N+2).*s(k+1)-b(l+1).*s(k+2));
  swap=s; s=t; t=swap;
end
j=floor(N/2):-1:0; s(j+2)=s(j+1);
for m=N-1:2*N-3,
  k=m+1-N:floor((m-1)/2); l=m-k; j=N-1-l;
  s(j+2)=cumsum(-(a(k+N+2)-a(l+1)).*t(j+2)-b(k+N+2).*s(j+2)+b(l+1).*s(j+3));
  j=j(length(j)); k=floor((m+1)/2);
  if rem(m,2)==0, a(k+N+2)=a(k+1)+(s(j+2)-b(k+N+2)*s(j+3))/t(j+3);
  else b(k+N+2)=s(j+2)/s(j+3);
  end
  swap=s; s=t; t=swap;
end
a(2*N+1)=a(N)-b(2*N+1)*s(2)/t(2);
ab=[a b];
*/
		static void r_kronrod(const IndexType N,VectorType & alpha, VectorType & beta)
		{
			assert(alpha.rows()==beta.rows());
			assert(alpha.rows() > 0);

			assert(alpha.rows() > ceil(3*N/2+1) );

			VectorType a=VectorType::Zero(2*N+1);
			VectorType b=VectorType::Zero(2*N+1);

			for(IndexType k=0;k<floor(3*N/2);++k)
			{
				a(k) = alpha(k)
			}

			for(IndexType k=0;k<ceil(3*N/2);++k)
			{
				b(k) =  beta(k)
			}

			VectorType s=VectorType::Zero(floor(N/2)+2);
			VectorType t=VectorType::Zero(floor(N/2)+2);

			t(1)=b(N+1) //t(2)=b(N+2)

			for(IndexType m=0;m<N-2;++m)
			{
				k=floor((m+1)/2):-1:0; l=m-k;
			}
		}

		
	};
}

#endif //LAURIE_HPP
