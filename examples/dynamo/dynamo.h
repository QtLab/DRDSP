#ifndef INCLUDED_DYNAMO
#define INCLUDED_DYNAMO
#include <DRDSP/dynamics/model.h>
#include <DRDSP/types.h>
#include <DRDSP/auto_diff.h>
#include <cmath>

namespace DRDSP {

	struct Dynamo : Model<> {

		double cAlpha = 1,
		       cOmega = -1000,
		       alphaB = 1;

		Dynamo() : Model<>( 2 * N ) {
			Init();
		}

		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,1> operator()( const MatrixBase<Derived>& x ) const {
			typedef Matrix<Derived::Scalar,-1,1> Col;
			auto a = compute_a( x );
			auto dsa = compute_ds( a );
			auto dta = compute_dt( a );
			auto ds2a = compute_ds2( a );
			auto dt2a = compute_dt2( a );
			auto b = compute_b( x );
			auto alpha = compute_alpha( a, b, dsa, dta );
			auto D2a = compute_D2( a, dsa, ds2a, dta, dt2a );
			auto aDot = compute_aDot( b, D2a, alpha );
			auto bDot = compute_bDot( a, b, dsa, dta, D2a, alpha );
			Col r(stateDim);
			for(uint32_t j=0;j<nJ;++j)
				r.segment(j*nI,nI) = aDot.col(j);
			for(uint32_t j=0;j<nJ;++j)
				r.segment(N+j*nI,nI) = bDot.col(j);
			return r;
		}
		
		SparseMatrix<double> Partials( const VectorXd& x ) const {
			return AutoDerivativeSparse( *this, x );
		}

		VectorXd InitialCondition() const;

		friend struct DynamoFamily;

	protected:
		static const uint32_t nI = 41;
		static const uint32_t nJ = 160;
		static const uint32_t N = nI*nJ;

		const double ds = 1.0 / double( nI - 1 ),
		             dth = (2.0 * PI) / double( nJ ),
		             eta0 = 1.0;

		VectorXd sinheta, cosheta, cotheta, sintheta, costheta, s;
		MatrixXd c, pi32, Bound1, Bound2, trig1;
		VectorXi jp1, jm1;

		void Init();
		
		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> compute_aDot( const MatrixBase<Derived>& b, const MatrixBase<Derived>& D2a, const MatrixBase<Derived>& alpha ) const {
			Matrix<Derived::Scalar,-1,-1> aDot = cAlpha * alpha.cwiseProduct( b ) + D2a;
			aDot.row(0).fill(aDot.row(1).mean());
			aDot.row(nI-1) = aDot.row(nI-2) * Bound1.transpose() + aDot.row(nI-3) * Bound2.transpose();
			return aDot;
		}
		
		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> compute_bDot( const MatrixBase<Derived>& a, const MatrixBase<Derived>& b, const MatrixBase<Derived>& dsa, const MatrixBase<Derived>& dta, const MatrixBase<Derived>& D2a, const MatrixBase<Derived>& alpha ) const {
			typedef Derived::Scalar Scalar;
			auto F = compute_F( dsa, dta );
			auto G = compute_G( a, dsa, dta, D2a, alpha );
			auto dsb = compute_ds( b );
			auto ds2b = compute_ds2( b );
			auto dtb = compute_dt( b );
			auto dt2b = compute_dt2( b );
			auto D2b = compute_D2( b, dsb, ds2b, dtb, dt2b );
			Matrix<Scalar,-1,-1> bDot = cOmega * F + cAlpha * G + D2b;
			bDot.row(0).fill(bDot.row(1).mean());
			bDot.row(nI-1).setZero();
			return bDot;
		}
		
		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> compute_a( const MatrixBase<Derived>& x ) const {
			Matrix<Derived::Scalar,-1,-1> a(nI,nJ);
			for(uint32_t j=0;j<nJ;++j)
				a.col(j) = x.segment(j*nI,nI);
			return a;
		}

		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> compute_b( const MatrixBase<Derived>& x ) const {
			Matrix<Derived::Scalar,-1,-1> b(nI,nJ);
			for(uint32_t j=0;j<nJ;++j)
				b.col(j) = x.segment(N+j*nI,nI);
			return b;
		}
		
		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> compute_alpha( const MatrixBase<Derived>& a, const MatrixBase<Derived>& b, const MatrixBase<Derived>& dsa, const MatrixBase<Derived>& dta ) const {
			typedef Derived::Scalar Scalar;
			auto mat1 = Matrix<Scalar,-1,-1>::Ones(nI,nJ);
			auto col1 = VectorXd::Ones(nI);
			auto normB2 = NormB2( a, b, dsa, dta );
			Matrix<Scalar,-1,-1> alpha = ( col1 * sintheta.transpose() ).cast<Scalar>().cwiseQuotient( mat1 + alphaB * normB2 );
			alpha.row(0).fill(alpha.row(1).mean());
			return alpha;
		}

		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> compute_ds( const MatrixBase<Derived>& x ) const {
			Matrix<Derived::Scalar,-1,-1> dsx;
			dsx.setZero(nI,nJ);
			for(uint32_t i=1;i<nI-1;++i)
				dsx.row(i) = x.row(i+1) - x.row(i-1);
			return dsx * ( 1.0 / (2.0*ds) );
		}
		
		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> compute_dt( const MatrixBase<Derived>& x ) const {
			Matrix<Derived::Scalar,-1,-1> dtx(nI,nJ);
			dtx.col(0) = x.col(1) - x.col(nJ-1);
			for(uint32_t j=1;j<nJ-1;++j)
				dtx.col(j) = x.col(j+1) - x.col(j-1);
			dtx.col(nJ-1) = x.col(0) - x.col(nJ-2);
			return dtx * ( 1.0 / (2.0*dth) );
		}

		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> compute_ds2( const MatrixBase<Derived>& x ) const {
			Matrix<Derived::Scalar,-1,-1> ds2x;
			ds2x.setZero(nI,nJ);
			for(uint32_t i=1;i<nI-1;++i)
				ds2x.row(i) = x.row(i+1) - 2.0 * x.row(i) + x.row(i-1);
			return ds2x * ( 1.0 / (ds*ds) );
		}
		
		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> compute_dt2( const MatrixBase<Derived>& x ) const {
			Matrix<Derived::Scalar,-1,-1> dt2x(nI,nJ);
			dt2x.col(0) = x.col(1) - 2.0 * x.col(0) + x.col(nJ-1);
			for(uint32_t j=1;j<nJ-1;++j)
				dt2x.col(j) = x.col(j+1) - 2.0 * x.col(j) + x.col(j-1);
			dt2x.col(nJ-1) = x.col(0) - 2.0 * x.col(nJ-1) + x.col(nJ-2);
			return dt2x * ( 1.0 / (dth*dth) );
		}

		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> compute_F( const MatrixBase<Derived>& dsa, const MatrixBase<Derived>& dta ) const {
			auto mat1 = MatrixXd::Ones(nI,nJ);
			return -1.5 * pi32.cwiseProduct(
				  dta.cwiseProduct( mat1 - cosheta * costheta.transpose() )
				- dsa.cwiseProduct( s.cwiseProduct(sinheta) * sintheta.transpose() )
			);
		}
		
		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> compute_G( const MatrixBase<Derived>& a, const MatrixBase<Derived>& dsa, const MatrixBase<Derived>& dta, const MatrixBase<Derived>& D2a, const MatrixBase<Derived>& alpha ) const {
			auto col1 = Matrix<double,-1,1>::Ones(nI);
			auto row1 = Matrix<double,1,-1>::Ones(nJ);
			auto dsalpha = compute_ds( alpha );
			auto dtalpha = compute_dt( alpha );
			auto cs = c.cwiseProduct( s * row1 );
			auto term1 = alpha.cwiseProduct( D2a );
			auto term2 = dsalpha.cwiseProduct( cs ).cwiseProduct(
						  a.cwiseProduct( trig1 )
						+ dsa.cwiseProduct( cs )
					);
			auto term3 = dtalpha.cwiseProduct( c ).cwiseProduct(
						  dta.cwiseProduct( c )
						- a.cwiseProduct( col1 * sintheta.transpose() )
					);
			
			return - term1 - term2 - term3;
		}
		
		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> Beta( const MatrixBase<Derived>& a, const MatrixBase<Derived>& dta ) const {
			typedef Derived::Scalar Scalar;
			auto col1 = VectorXd::Ones(nI);
			return c.cwiseProduct(dta) - a.cwiseProduct( col1 * sintheta.transpose() );
		}

		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> Btheta( const MatrixBase<Derived>& a, const MatrixBase<Derived>& dsa ) const {
			auto row1 = Matrix<double,1,-1>::Ones(nJ);
			return -a.cwiseProduct(trig1) + c.cwiseProduct(dsa).cwiseProduct( s * row1 );
		}

		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> NormB2( const MatrixBase<Derived>& a, const MatrixBase<Derived>& b, const MatrixBase<Derived>& dsa, const MatrixBase<Derived>& dta ) const {
			auto b1 = Beta(a,dta);
			auto b2 = Btheta(a,dsa);
			return b1.cwiseProduct(b1) + b2.cwiseProduct(b2) + b.cwiseProduct(b);
		}
		
		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,-1> compute_D2( const MatrixBase<Derived>& x, const MatrixBase<Derived>& dsx, const MatrixBase<Derived>& ds2x, const MatrixBase<Derived>& dtx, const MatrixBase<Derived>& dt2x ) const {
			auto c2 = c.array().square().matrix().eval();
			auto s2 = s.array().square().matrix().eval();
			auto col1 = Matrix<double,-1,1>::Ones(nI);
			auto row1 = Matrix<double,1,-1>::Ones(nJ);
			auto term1 = ds2x.cwiseProduct( s2 * row1 );
			auto term2 = dsx.cwiseProduct( ( s * row1 ).cwiseProduct( ( col1 - cotheta ) * row1 + ( sinheta * row1 ).cwiseQuotient( c ) ) );
			auto term4 = dtx.cwiseProduct( ( col1 * sintheta.transpose() ).cwiseQuotient( c ) );
			auto term5 = x.cwiseQuotient( sinheta.array().square().matrix() * row1 );
			return c2.cwiseProduct(
				  term1
				+ term2
				+ dt2x
				- term4
				- term5
			);
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

