#ifndef INCLUDED_DYNAMICS_RADIAL_BASIS
#define INCLUDED_DYNAMICS_RADIAL_BASIS
#include "../types.h"
#include "../auto_diff.h"
#include "polyharmonic_spline.h"

namespace DRDSP {
	
	using std::sqrt;
	using std::exp;

	struct Gaussian  {
		double scale = 1.0;

		Gaussian() = default;

		explicit Gaussian( double scale ) : scale(scale) {}

		template<typename T>
		T operator()( T r ) const {
			T x = scale * r;
			return exp(-x*x);
		}

		template<typename T>
		T Derivative( T r ) const {
			T x = scale * r;
			return (-2.0 * scale) * x * exp(-x*x);
		}
	};

	struct Multiquadratic {
		double scale = 1.0;
		
		Multiquadratic() = default;

		explicit Multiquadratic( double scale ) : scale(scale) {}

		template<typename T>
		T operator()( T r ) const {
			T x = scale * r;
			return sqrt( T(1) + x*x );
		}

		template<typename T>
		T Derivative( T r ) const {
			T x = scale * r;
			return ( scale * x ) / sqrt( T(1) + x*x );
		}
	};

	struct InverseQuadratic {
		double scale = 1.0;
		
		InverseQuadratic() = default;
		
		explicit InverseQuadratic( double scale ) : scale(scale) {}
		
		template<typename T>
		T operator()( T r ) const {
			T x = scale * r;
			return T(1) / ( T(1) + x*x );
		}

		template<typename T>
		T Derivative( T r ) const {
			T x = scale * r;
			T y = T(1) + x*x;
			return ( (-2.0 * scale) * x ) / ( y*y );
		}
	};

	struct InverseMultiquadratic {
		double scale = 1.0;

		InverseMultiquadratic() = default;

		explicit InverseMultiquadratic( double scale ) : scale(scale) {}

		template<typename T>
		T operator()( T r ) const {
			T x = scale * r;
			return T(1) / sqrt( T(1) + x*x );
		}

		template<typename T>
		T Derivative( T r ) const {
			T x = scale * r;
			T y = T(1) + x*x;
			return ( -scale * x ) / sqrt( y*y*y );
		}
	};

	template<typename F>
	struct RBF {
		typedef F RadialType;
		VectorXd weight, centre;
		F func;
	
		RBF() = default;
		
		explicit RBF( const F& f ) : func(f) {}

		VectorXd operator()( const VectorXd& x ) const {
			return weight * func( (x-centre).norm() );
		}

		MatrixXd Derivative( const VectorXd& x ) const {
			VectorXd r = x - centre;
			double rnorm = r.norm();
			if( rnorm == 0.0 ) return VectorXd::Zero(x.size());
			return weight * ( ( func.Derivative( rnorm ) / rnorm ) * r ).transpose();
		}

		MatrixXd LinearWeight( const VectorXd& x ) const {
			return MatrixXd::Identity(x.size(),x.size()) * func( (x-centre).norm() );
		}
		
		MatrixXd LinearDerivativeWeight( const VectorXd& x ) const {
			int64_t dim = x.size();
			VectorXd r = x - centre;
			double rnorm = r.norm();
			if( rnorm == 0.0 ) return MatrixXd::Zero( dim * dim, dim );
			MatrixXd C( dim * dim, dim );
			for(int64_t i=0;i<dim;++i) {
				C.block(i*dim,0,dim,dim).setIdentity() *= r[i];
			}
			return C * ( func.Derivative( rnorm ) / rnorm );
		}
	};

	template<typename F>
	struct EquiRBFZ2 {
		typedef F RadialType;
		VectorXd weight, centre;
		F func;
	
		EquiRBFZ2() = default;
		
		explicit EquiRBFZ2( const F& f ) : func(f) {}
		
		VectorXd operator()( const VectorXd& x ) const {
			return weight * ( func( (x-centre).norm() ) - func( (x+centre).norm() ) );
		}

		MatrixXd Derivative( const VectorXd& x ) const {
			VectorXd sum;
			sum.setZero( x.size() );
			VectorXd r = x - centre;
			double rnorm = r.norm();
			if( rnorm != 0.0 ) {
				sum += ( func.Derivative( rnorm ) / rnorm ) * r;
			}
			r = x + centre;
			rnorm = r.norm();
			if( rnorm != 0.0 ) {
				sum -= ( func.Derivative( rnorm ) / rnorm ) * r;
			}
			return weight * sum.transpose();
		}

		MatrixXd LinearWeight( const VectorXd& x ) const {
			return MatrixXd::Identity(x.size(),x.size()) * ( func( (x-centre).norm() ) - func( (x+centre).norm() ) );
		}

		MatrixXd LinearDerivativeWeight( const VectorXd& x ) const {
			int64_t dim = x.size();
			MatrixXd C( dim * dim, dim );
			C.setZero();
			VectorXd r = x - centre;
			double rnorm = r.norm();
			if( rnorm != 0.0 ) {
				double dphi = func.Derivative( rnorm ) / rnorm;
				for(int64_t i=0;i<dim;++i) {
					C.block(i*dim,0,dim,dim) += MatrixXd::Identity(dim,dim) * ( r[i] * dphi );
				}
			}
			r = x + centre;
			rnorm = r.norm();
			if( rnorm != 0.0 ) {
				double dphi = func.Derivative( rnorm ) / rnorm;
				for(int64_t i=0;i<dim;++i) {
					C.block(i*dim,0,dim,dim) -= MatrixXd::Identity(dim,dim) * ( r[i] * dphi );
				}
			}
			return C;
		}
	};

	template<typename F,int N>
	struct EquiRBFCyclic {
		typedef F RadialType;
		VectorXd weight, centre;
		MatrixXd generator;
		F func;

		EquiRBFCyclic() = default;
		
		explicit EquiRBFCyclic( const F& f ) : func(f) {}
		
		VectorXd operator()( const VectorXd& x ) const {
			VectorXd sum;
			sum.setZero( x.size() );
			VectorXd c = centre, w = weight;
			for(int i=0;i<N;++i) {
				w = generator * w;
				c = generator * c;
				sum.noalias() += w * func( (x-c).norm() );
			}
			return sum;
		}

		MatrixXd Derivative( const VectorXd& x ) const {
			VectorXd r;
			MatrixXd sum;
			sum.setZero( x.size() );
			VectorXd c = centre, w = weight;
			for(int i=0;i<N;++i) {
				w = generator * w;
				c = generator * c;
				r = x - c;
				double rnorm = r.norm();
				if( rnorm == 0.0 ) continue;
				sum += w * ( ( func.Derivative( rnorm ) / rnorm ) * r ).transpose();
			}
			return sum;
		}

		MatrixXd LinearWeight( const VectorXd& x ) const {
			MatrixXd sum( generator.rows(), generator.cols() );
			MatrixXd g = MatrixXd::Identity( generator.rows(), generator.cols() );
			sum.setIdentity() *= func( (x-c).norm() );
			for(int i=1;i<N;++i) {
				g = generator * g;
				sum.noalias() += g * func( (x-g*c).norm() );
			}
			return sum;
		}

		MatrixXd LinearDerivativeWeight( const VectorXd& x ) const {
			int64_t dim = x.size();
			MatrixXd C( dim * dim, dim );
			C.setZero();
			VectorXd r = x - centre;
			double rnorm = r.norm();
			if( rnorm != 0.0 ) {
				double dphi = func.Derivative( rnorm ) / rnorm;
				for(int64_t i=0;i<dim;++i) {
					C.block(i*dim,0,dim,dim).setIdentity() *= dphi * r[i];
				}
			}
			MatrixXd g = MatrixXd::Identity( generator.rows(), generator.cols() );
			for(int i=1;i<N;++i) {
				g = generator * g;
				r = x - g * centre;
				rnorm = r.norm();
				if( rnorm == 0.0 ) continue;
				double dphi = func.Derivative( rnorm ) / rnorm;
				for(int64_t i=0;i<dim;++i) {
					C.block(i*dim,0,dim,dim) += g * ( dphi * r[i] );
				}
			}
			return C;
		}
	};

	template<typename F>
	struct EquiRBFFinite {
		typedef F RadialType;
		VectorXd weight, centre;
		vector<MatrixXd> group;
		F func;

		EquiRBFFinite() = default;
		
		explicit EquiRBFFinite( const F& f ) : func(f) {}

		VectorXd operator()( const VectorXd& x ) const {
			VectorXd sum;
			sum.setZero( x.size() );
			for( const auto& g : group ) {
				sum += g * weight * func( ( x - g * centre ).norm() );
			}
			return sum;
		}

		MatrixXd Derivative( const VectorXd& x ) const {
			VectorXd r;
			MatrixXd sum;
			sum.setZero( x.size() );
			for( const auto& g : group ) {
				r = x - g * centre;
				double rnorm = r.norm();
				if( rnorm == 0.0 ) continue;
				sum.noalias() += g * weight * ( ( func.Derivative( rnorm ) / rnorm ) * r ).transpose();
			}
			return sum;
		}

		MatrixXd LinearWeight( const VectorXd& x ) const {
			MatrixXd sum;
			sum.setZero( generator.rows(), generator.cols() );
			for( const auto& g : group ) {
				sum += g * func( ( x - g * centre ).norm() );
			}
			return sum;
		}

		MatrixXd LinearDerivativeWeight( const VectorXd& x ) const {
			int64_t dim = x.size();
			MatrixXd C;
			C.setZero( dim * dim, dim );
			VectorXd r;
			for( const auto& g : group ) {
				r = x - g * centre;
				double rnorm = r.norm();
				if( rnorm == 0.0 ) continue;
				double dphi = func.Derivative( rnorm ) / rnorm;
				for(int64_t i=0;i<dim;++i) {
					C.block(i*dim,0,dim,dim) += g * ( dphi * r[i] );
				}
			}
			return C;
		}
	};

	template<typename F>
	struct EquiRBFSO2;

	template<>
	struct EquiRBFSO2<InverseQuadratic> {
		typedef InverseQuadratic RadialType;
		VectorXd weight, centre;
		Matrix<double,-1,2> W;
		InverseQuadratic func;
	
		EquiRBFSO2() = default;
		
		explicit EquiRBFSO2( const InverseQuadratic& f ) : func(f) {}

		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,1> operator()( const MatrixBase<Derived>& x ) const {
			typedef Derived::Scalar Scalar;
			Matrix<Scalar,2,1> Wtx = W.adjoint() * x;
			Matrix<double,2,1> Wtc = W.adjoint() * centre;
			Scalar rR = Wtx.squaredNorm() * Wtc.squaredNorm();
			if( rR == Scalar() ) {
				return ( 2.0 * PI * func( (x-centre.cast<Scalar>()).norm() ) ) * ( weight - W * ( W.adjoint() * weight ) ).cast<Scalar>();
			}
			auto rotation = ComputeRotation( Wtc, Wtx ).transpose();
			return (-PI / sqrt(rR)) * ( W.cast<Scalar>() * ( rotation * ( W.adjoint() * weight ) ) );
		}

		MatrixXd Derivative( const VectorXd& x ) const {
			return AutoDerivative( *this, x );
		}

		MatrixXd LinearWeight( const VectorXd& x ) const {
			Vector2d Wtx = W.adjoint() * x;
			Vector2d Wtc = W.adjoint() * centre;
			double rR = Wtx.squaredNorm() * Wtc.squaredNorm();
			if( rR == 0.0 ) {
				return ( 2.0 * PI * func( (x-centre).norm() ) ) * ( MatrixXd::Identity(W.rows(),W.rows()) - W * W.adjoint() );
			}
			auto rotation = ComputeRotation( Wtc, Wtx ).transpose();
			return (-PI / sqrt(rR)) * ( W * rotation * W.adjoint() );
		}

		MatrixXd LinearDerivativeWeight( const VectorXd& x ) const {
			int64_t dim = x.size();
			MatrixXd C;
			C.setZero( dim * dim, dim );
			Vector2d Wtx = W.adjoint() * x;
			Vector2d Wtc = W.adjoint() * centre;
			double rR = Wtx.squaredNorm() * Wtc.squaredNorm();
			if( rR == 0.0 ) {
				VectorXd r = x - centre;
				double rnorm = r.norm();
				if( rnorm == 0.0 ) return C;
				double dphi = 2.0 * PI * func.Derivative( rnorm ) / rnorm;
				MatrixXd g = MatrixXd::Identity(dim,dim) - W * W.adjoint();
				for(int64_t i=0;i<dim;++i) {
					C.block(i*dim,0,dim,dim) = g * ( dphi * r[i] );
				}
				return C;
			}

			Matrix2d rotation = ComputeRotation( Wtc, Wtx );
			
			// todo

			return C;
		}

	protected:

		template<typename Derived1,typename Derived2>
		static Matrix<typename internal::scalar_product_traits<typename Derived1::Scalar,typename Derived2::Scalar>::ReturnType,2,2> ComputeRotation( const MatrixBase<Derived1>& x, const MatrixBase<Derived2>& y ) {
			typedef typename internal::scalar_product_traits<typename Derived1::Scalar,typename Derived2::Scalar>::ReturnType Scalar;
			Scalar costheta = x.dot(y) / ( x.norm() * y.norm() );
			Scalar sintheta = sqrt( Scalar(1) - costheta*costheta );
			Matrix<Scalar,2,2> rotation;
			rotation(0,0) = costheta;
			rotation(1,0) = sintheta;
			rotation(0,1) = -sintheta;
			rotation(1,1) = costheta;
			return rotation;
		}
	};
	
}

#endif
