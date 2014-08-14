#ifndef LAURIE_HPP
#define LAURIE_HPP

namespace gauss_konrad {

	template<typename _RealType>
	class Laurie
	{
	public:
		typedef _RealType RealType;
		typedef typename Eigen::Matrix<RealType,Eigen::Dynamic,1>  VectorType;
		typedef typename Eigen::Matrix<RealType,Eigen::Dynamic,Eigen::Dynamic> MatrixType;
		typedef typename VectorType::Index IndexType;
		typedef typename Eigen::SelfAdjointEigenSolver<MatrixType> SelfAdjointEigenSolverType;

		static void r_jacobi(const IndexType N,const RealType a,const RealType b,VectorType & a_out, VectorType & b_out)
		{
			assert(a > -1);
			assert(b > -1);

			assert(a_out.rows()==b_out.rows());
			assert(a_out.rows() > 0);

			assert(N<=a_out.rows());

			a_out(0) = (b-a)/(a+b+2.);
			b_out(0) = pow(2.,(a+b+1.))*gamma(a+1.)*gamma(b+1.)/gamma(a+b+2.);

			for(IndexType n=1;n<N;++n)
			{
				RealType nab = 2.*n+a+b;
				a_out(n) = (b*b - a*a)/(nab*(nab+2.));
				b_out(n) =  4.*(n+a)*(n+b)*n*(n+a+b)/(nab*nab*(nab+1.)*(nab-1.));
			}
		}

		static void r_jacobi_01(const IndexType N,const RealType a,const RealType b,VectorType & a_out, VectorType & b_out)
		{
			assert(a > -1);
			assert(b > -1);

			assert(a_out.rows()==b_out.rows());
			assert(a_out.rows() > 0);

			assert(N<=a_out.rows());

			r_jacobi(N,a, b, a_out, b_out);

			for(IndexType n=0;n<N;++n)
			{
				a_out(n)  = (1+a_out(n))/2.;
			}

			b_out(0) = b_out(0)/pow(2,a+b+1.);

			for(IndexType n=1;n<N;++n)
			{
				b_out(n) = b_out(n)/4.;
			}

		}


		static void r_kronrod(const IndexType N,VectorType const & a_in, VectorType const & b_in,
			VectorType & a, VectorType & b)
		{
			assert(a_in.rows()==b_in.rows());
			assert(a_in.rows()>0);
			assert(a_in.rows() >= ceil(3*N/2)+1 );

			assert(a.rows()==2*N+1);
			assert(b.rows()==2*N+1);

			a=VectorType::Zero(2*N+1);
			b=VectorType::Zero(2*N+1);

			for(IndexType k=0;k<floor(3*N/2)+1;++k)
			{
				a(k) = a_in(k);
			}

			for(IndexType k=0;k<ceil(3*N/2)+1;++k)
			{
				b(k) = b_in(k);
			}

			VectorType s=VectorType::Zero(floor(N/2)+2);
			VectorType t=VectorType::Zero(floor(N/2)+2);

			t(1)=b(N+1);

			

			for(IndexType m=0;m<N-2+1;++m)
			{
				RealType u=0;
				for(IndexType k=floor((m+1)/2);k>=0;--k)
				{
					IndexType l=m-k;
					u = u + ( a(k+N+1)-a(l) )*t(k+1) + b(k+N+1)*s(k) - b(l)*s(k+1);
					s(k+1) = u;
				}
				VectorType swap=s;
				s=t;
				t=swap;
			}

			for(IndexType j=floor(N/2);j>=0;--j)
			{
				s(j+1)=s(j);
			}

			for(IndexType m=N-1;m<2*N-3+1;++m)
			{
				IndexType k;
				IndexType j;
				RealType u=0;
				for(k=m+1-N;k<floor((m-1)/2)+1;++k)
				{
					IndexType l=m-k;
					j=N-1-l;
					//std::cout<<j<<"\t";
					//s(j+2)=cumsum(-(a(k+N+2)-a(l+1)).*t(j+2)-b(k+N+2).*s(j+2)+b(l+1).*s(j+3));
					u = u - ( a(k+N+1)-a(l) )*t(j+1) - b(k+N+1)*s(j+1) + b(l)*s(j+2);
					s(j+1) = u;
				}

				//std::cout<<std::endl;

				//std::cout<<s<<"\n"<<std::endl;

				k=floor((m+1)/2);

				std::cout<<k<<std::endl;

				if(m % 2 == 0)
				{
					//a(k+N+2)=a(k+1)+(s(j+2)-b(k+N+2)*s(j+3))/t(j+3);
					a(k+N+1)=a(k)+(s(j+1)-b(k+N+1)*s(j+2))/t(j+2);
					std::cout<<"a(k+N+1) 2 "<<a(k+N+1)<<std::endl;
				}
				else
				{
					//b(k+N+2)=s(j+2)/s(j+3);
					b(k+N+1)=s(j+1)/s(j+2);
				}
			}

			a(2*N)=a(N-1)-b(2*N)*s(1)/t(1);
		}

		/*
		static void r_kronrod(const IndexType N,VectorType & a, VectorType & b)
		{
			assert(a.rows()==b.rows());
			assert(a.rows() > 0);
			assert(N>0);
			assert(a.rows() == 2*N+1 );

			VectorType s=VectorType::Zero(floor(N/2)+2);
			VectorType t=VectorType::Zero(floor(N/2)+2);

			t(1)=b(N+1);

			for(IndexType m=1;m<=N-2+1;++m)
			{
				RealType u=0;
				//std::cout<<"m= "<<m<<std::endl;
				for(IndexType k=floor((m+1)/2);k>0;--k)
				{
					IndexType l=m-k;
					u += ( a(k+N+1) - a(l) )*t(k) + b(k+N+1)*s(k-1) - b(l)*s(k);
					s(k) = u;
				}

				VectorType swap(s);
				s=t;
				t=swap;
			}

			return;

			for(IndexType j=floor(N/2);j>0;--j)
			{
				s(j) = s(j-1);
			}

			return;

			for(IndexType m=N-1;m<=2*N-3;++m)
			{
				RealType u=0;
				IndexType j,k;
				for(k=m+1-N;k<=floor((m-1)/2);++k)
				{
					k += 1;
					IndexType l = m-k;
					j = N-1-l;
					u -= ( a(k+N+1) -  a(l) )*t(j) - b(k+N+1)*s(j) + b(l)*s(j+1);
					s(j) = u;
				}

				if(m % 2 == 0)
				{
					k = m/2;
					a(k+N+1) = a(k) + ( s(j) -b(k+N+1)*s(j+1) )/t(j+1);
				}
				else
				{
					k = (m+1)/2;
					b(k+N+1) =  s(j)/s(j+1);
				}

				VectorType swap(s);
				s=t;
				t=swap;
			}

			a(2*N) = a(N-1) - b(2*N)*s(1)/t(1);
		}
		*/
		/*
		static void kronrod(const IndexType N,VectorType const & alpha, VectorType const & beta,VectorType & x, VectorType & w)
		{
			//CHECK NEEDED LIKE THE ONE ON LINE 21 IN KONRAD.M

			assert(alpha.rows()==2*N+1)
			assert(alpha.rows()==beta.rows());

			MatrixType J=MatrixType::Zero(2*N+1,2*N+1);

			for(IndexType k=0;k<2*N;++k)
			{
				J(k,k)=alpha(k);
				J(k,k+1)=sqrt(beta(k+1));
				J(k+1,k)=J(k,k+1);
			}

			J(2*N,2*N)=alpha(2*N);

			SelfAdjointEigenSolverType es(J);

			assert(es.info()==Eigen::Success);

			VectorType x=es.eigenvalues();
			MatrixType V=es.eigenvectors();

			w=beta(0)*(V.row(1).array()*V.row(1).array()).matrix();

			//the elements for d are ordered in the increasing order
			//so we don't need to sort this?
		}

		static void mpkonrad(const IndexType N,VectorType & x, VectorType & w)
		{
			VectorType alpha(2*N+1);
			VectorType beta(2*N+1);
			r_jacobi_01( 2*N, RealType(0), RealType(0), alpha, beta);

			VectorType x(2*N+1);
			VectorType w(2*N+1);

			kronrod(N,alpha,beta,x,w);

			x = 2.*x-1.;
			w = 2.*w;

		}
		*/
		
	};
}

#endif //LAURIE_HPP
