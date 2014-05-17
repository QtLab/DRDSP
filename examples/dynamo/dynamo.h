#ifndef INCLUDED_DYNAMO
#define INCLUDED_DYNAMO
#include <DRDSP/dynamics/model.h>
#include <DRDSP/types.h>
#include <DRDSP/auto_diff.h>
#include <cmath>

namespace DRDSP {

	struct Dynamo : Model<> {

		static const uint32_t nI = 41;
		static const uint32_t nJ = 160;
		static const uint32_t N = nI*nJ;
		double cAlpha, cOmega, alphaB;

		Dynamo() :
			Model<>( 2 * N ),
			ds( 1.0 / double( nI - 1 ) ),
			dth( (2.0 * M_PI) / double( nJ ) ),
			cAlpha( 1 ),
			cOmega( -1000 ),
			alphaB( 1.0 ),
			eta0( 1.0 )
		{
			Init();
		}

		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,1> operator()( const MatrixBase<Derived>& x ) const {

			Matrix<typename Derived::Scalar,-1,-1> alphaPC, aDotPC, bDotPC;
			
			alphaPC.setZero(nI,nJ);
			for(uint32_t j=0;j<nJ;j++)
				for(uint32_t i=1;i<nI-1;i++)
					alphaPC(i,j) = sintheta(j) / ( 1.0 + alphaB * NormB2(x,i,j) );
		
			Derived::Scalar t(0);
			for(uint32_t k=0;k<nJ;k++)
				t += alphaPC(1,k);
			t /= nJ;
			for(uint32_t j=0;j<nJ;j++)
				alphaPC(0,j) = t;

			aDotPC.setZero(nI,nJ);
			for(uint32_t j=0;j<nJ;j++)
				for(uint32_t i=1;i<nI-1;i++)
					aDotPC(i,j) = cAlpha * alphaPC(i,j) * b(x,i,j) + D2a(x,i,j);

			for(uint32_t j=0;j<nJ;j++)
				for(uint32_t k=0;k<nJ;k++)
					aDotPC(nI-1,j) += Bound1(j,k) * aDotPC(nI-2,k) + Bound2(j,k) * aDotPC(nI-3,k);
		
			t = 0.0;
			for(uint32_t k=0;k<nJ;k++)
				t += aDotPC(1,k);
			t /= nJ;
			for(uint32_t j=0;j<nJ;j++)
				aDotPC(0,j) = t;

			bDotPC.setZero(nI,nJ);
			for(uint32_t j=0;j<nJ;j++)
				for(uint32_t i=1;i<nI-1;i++)
					bDotPC(i,j) = cOmega * F(x,i,j) + cAlpha * G(x,alphaPC,i,j) + D2b(x,i,j);
		
			t = 0.0;
			for(uint32_t k=0;k<nJ;k++)
				t += bDotPC(1,k);
			t /= nJ;
			for(uint32_t j=0;j<nJ;j++)
				bDotPC(0,j) = t;

			Matrix<typename Derived::Scalar,-1,1> r(dimension);
			uint32_t k=0;
			for(uint32_t i=0;i<nI;i++) {
				for(uint32_t j=0;j<nJ;j++)
					r(k++) = aDotPC(i,j);
			}
			for(uint32_t i=0;i<nI;i++) {
				for(uint32_t j=0;j<nJ;j++)
					r(k++) = bDotPC(i,j);
			}
			return r;
		}

		SparseMatrix<double> Partials( const VectorXd& x ) const {
			return AutoDerivativeSparse( *this, x );
		}

		VectorXd InitialCondition() const;

	protected:
		const double ds, dth, eta0;
		VectorXd sinheta, cosheta, cotheta, sintheta, costheta;
		MatrixXd c, pi32, Bound1, Bound2, trig1;
		VectorXi jp1, jm1;

		void Init();

		template<typename Derived>
		typename Derived::Scalar a( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return x(nJ*i+j);
		}

		template<typename Derived>
		typename Derived::Scalar b( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return x(N+nJ*i+j);
		}

		template<typename Derived>
		typename Derived::Scalar dsa( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return ( a(x,i+1,j) - a(x,i-1,j) ) / ( 2.0*ds );
		}

		template<typename Derived>
		typename Derived::Scalar dsb( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return ( b(x,i+1,j) - b(x,i-1,j) ) / ( 2.0*ds );
		}

		template<typename Derived>
		typename Derived::Scalar dta( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return ( a(x,i,jp1(j)) - a(x,i,jm1(j)) ) / ( 2.0*dth );
		}

		template<typename Derived>
		typename Derived::Scalar dtb( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return ( b(x,i,jp1(j)) - b(x,i,jm1(j)) ) / ( 2.0*dth );
		}

		//template<typename Derived>
		//typename Derived::Scalar dsta( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
		//	return ( a(x,i+1,jp1(j)) - a(x,i+1,jm1(j)) - a(x,i-1,jp1(j)) + a(x,i-1,jm1(j)) ) / ( 4.0*ds*dth );
		//}

		//template<typename Derived>
		//typename Derived::Scalar dstb( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
		//	return ( b(x,i+1,jp1(j)) - b(x,i+1,jm1(j)) - b(x,i-1,jp1(j)) + b(x,i-1,jm1(j)) ) / ( 4.0*ds*dth );
		//}

		template<typename Derived>
		typename Derived::Scalar ds2a( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return ( a(x,i+1,j) - 2.0*a(x,i,j) + a(x,i-1,j) ) / ( ds*ds );
		}

		template<typename Derived>
		typename Derived::Scalar ds2b( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return ( b(x,i+1,j) - 2.0*b(x,i,j) + b(x,i-1,j) ) / ( ds*ds );
		}

		template<typename Derived>
		typename Derived::Scalar dt2a( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return ( a(x,i,jp1(j)) - 2.0*a(x,i,j) + a(x,i,jm1(j)) ) / ( dth*dth );
		}

		template<typename Derived>
		typename Derived::Scalar dt2b( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return ( b(x,i,jp1(j)) - 2.0*b(x,i,j) + b(x,i,jm1(j)) ) / ( dth*dth );
		}
		
		double theta( uint32_t j ) const {
			return dth * j;
		}

		double s( uint32_t i ) const {
			return ds * i;
		}

		double eta( uint32_t i ) const {
			return eta0 - ::log(s(i));
		}
		/*
		double c( uint32_t i, uint32_t j ) const {
			return cosheta(i) - costheta(j);
		}

		double pi( uint32_t i, uint32_t j ) const {
			return sinheta(i) / c(i,j);
		}*/

		template<typename Derived>
		typename Derived::Scalar F( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return -1.5 * pi32(i,j) * ( (1.0-cosheta(i)*costheta(j)) * dta(x,i,j) - s(i) * sinheta(i) * sintheta(j) * dsa(x,i,j) );
		}
		
		template<typename Derived,typename Derived2>
		typename Derived::Scalar G( const MatrixBase<Derived>& x, const MatrixBase<Derived2>& alpha, uint32_t i, uint32_t j ) const {
			return -alpha(i,j)*D2a(x,i,j) - c(i,j)*s(i) * dsalphaApprox(alpha,i,j) * ( trig1(i,j)*a(x,i,j) + c(i,j)*s(i)*dsa(x,i,j) ) - c(i,j)*dtalphaApprox(alpha,i,j)*( c(i,j)*dta(x,i,j) - a(x,i,j)*sintheta(j) );
		}
		
		template<typename Derived>
		typename Derived::Scalar Beta( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return c(i,j)*dta(x,i,j) - a(x,i,j)*sintheta(j);
		}

		template<typename Derived>
		typename Derived::Scalar Btheta( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return -a(x,i,j)*trig1(i,j) + c(i,j)*s(i)*dsa(x,i,j);
		}

		template<typename Derived>
		typename Derived::Scalar Bphi( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return b(x,i,j);
		}

		template<typename Derived>
		typename Derived::Scalar NormB2( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			Derived::Scalar b1 = Beta(x,i,j);
			Derived::Scalar b2 = Btheta(x,i,j);
			Derived::Scalar b3 = Bphi(x,i,j);
			return b1*b1 + b2*b2 + b3*b3;
		}
		
		template<typename Derived>
		typename Derived::Scalar D2a( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return c(i,j)*c(i,j) * ( s(i)*s(i) * ds2a(x,i,j) + s(i) * (1.0 - cotheta(i) + sinheta(i)/c(i,j)) * dsa(x,i,j) + dt2a(x,i,j) - (sintheta(j)/c(i,j)) * dta(x,i,j) - a(x,i,j)/(sinheta(i)*sinheta(i)) );
		}
		
		template<typename Derived>
		typename Derived::Scalar D2b( const MatrixBase<Derived>& x, uint32_t i, uint32_t j ) const {
			return c(i,j)*c(i,j) * ( s(i)*s(i) * ds2b(x,i,j) + s(i) * (1.0 - cotheta(i) + sinheta(i)/c(i,j)) * dsb(x,i,j) + dt2b(x,i,j) - (sintheta(j)/c(i,j)) * dtb(x,i,j) - b(x,i,j)/(sinheta(i)*sinheta(i)) );
		}

		template<typename Derived>
		typename Derived::Scalar dsalphaApprox( const MatrixBase<Derived>& alpha, uint32_t i, uint32_t j ) const {
			return ( alpha(i+1,j) - alpha(i-1,j) ) / (2.0*ds);
		}

		template<typename Derived>
		typename Derived::Scalar dtalphaApprox( const MatrixBase<Derived>& alpha, uint32_t i, uint32_t j ) const {
			return ( alpha(i,jp1(j)) - alpha(i,jm1(j)) ) / (2.0*dth);
		}

	};

	struct DynamoFamily : Family<Dynamo> {

		DynamoFamily() : Family<Dynamo>(2*Dynamo::N,1) {}

		Dynamo operator()( const VectorXd& parameter ) const {
			Dynamo model;
			model.cAlpha = parameter[0];
			return model;
		}
	
	};

}

#endif

