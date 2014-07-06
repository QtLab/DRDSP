#ifndef INCLUDED_DYNAMICS_RADIAL_BASIS
#define INCLUDED_DYNAMICS_RADIAL_BASIS
#include "../types.h"
#include "../auto_diff.h"
#include <cmath>

namespace DRDSP {
	
	using std::sqrt;
	using std::exp;
	using std::log;

	struct ThinPlateSpline {
		template<typename T>
		T operator()( T r ) const {
			return r*r*log(r);
		}
		template<typename T>
		T Derivative( T r ) const {
			return r * ( T(1) + 2.0 * log(r) );
		}
	};

	struct PolyharmonicSpline3 {
		template<typename T>
		T operator()( T r ) const {
			return r*r*r;
		}
		template<typename T>
		T Derivative( T r ) const {
			return T(3)*r*r;
		}
	};

	struct Gaussian  {
		double scale;

		Gaussian() : Gaussian(1.0) {}

		explicit Gaussian( double scale ) : scale(scale) {}

		template<typename T>
		T operator()( T r ) const {
			T x = scale * r;
			return exp(-x*x);
		}

		template<typename T>
		T Derivative( T r ) const {
			T x = scale * r;
			return -T(2) * scale * x * exp(-x*x);
		}
	};

	struct Multiquadratic {
		double scale;
		
		Multiquadratic() : Multiquadratic(1.0) {}

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
		double scale;
		
		InverseQuadratic() : InverseQuadratic(1.0) {}
		
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
			return ( -2.0 * scale * x ) / ( y*y );
		}
	};

	struct InverseMultiquadratic {
		double scale;
		
		InverseMultiquadratic() : InverseMultiquadratic(1.0) {}

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
		VectorXd centre;
		F func;
	
		RBF() = default;
		
		explicit RBF( const F& f ) : func(f) {}
		
		explicit RBF( const VectorXd& centre ) : centre(centre) {}
		
		RBF( const F& f, const VectorXd& centre ) : func(f), centre(centre) {}

		double operator()( const VectorXd& x ) const {
			return func( (x-centre).norm() );
		}

		VectorXd Derivative( const VectorXd& x ) const {
			VectorXd r = x - centre;
			double rnorm = r.norm();
			if( rnorm == 0.0 ) return VectorXd::Zero(x.size());
			return ( func.Derivative( rnorm ) / rnorm ) * r;
		}
	};

	template<typename F,typename G>
	struct EquiRBFVecField {
		VectorXd weight, centre;
		F func;
		G group;
	
		EquiRBFVecField() = default;
		
		explicit EquiRBFVecField( const F& f ) : func(f) {}
		
		explicit EquiRBFVecField( const VectorXd& centre ) : centre(centre) {}

		EquiRBFVecField( const F& f, const G& g ) : func(f), group(g) {}

		EquiRBFVecField( const F& f, const VectorXd& centre ) : func(f), centre(centre) {}
		
		EquiRBFVecField( const F& f, const VectorXd& weight, const VectorXd& centre ) : func(f), weight(weight), centre(centre) {}

		EquiRBFVecField( const F& f, const VectorXd& weight, const VectorXd& centre, const G& g ) : func(f), weight(weight), centre(centre), group(g) {}

		VectorXd operator()( const VectorXd& x ) const {
			VectorXd sum;
			sum.setZero( x.size() );
			for( const auto& g : group ) {
				sum += g * weight * func( ( x - g * centre ).norm() );
			}
			return sum;
		}

		MatrixXd Partials( const VectorXd& x ) const {
			VectorXd r;
			MatrixXd sum;
			sum.setZero( x.size() );
			for( const auto& g : group ) {
				r = x - g * centre;
				double rnorm = r.norm();
				if( rnorm == 0.0 ) continue;
				sum += g * weight * ( ( func.Derivative( rnorm ) / rnorm ) * r ).transpose();
			}
			return sum;
		}
	};

	template<typename F>
	struct EquiRBFSO2;

	template<>
	struct EquiRBFSO2<InverseQuadratic> {
		VectorXd weight, centre;
		Matrix<double,-1,2> W;
		InverseQuadratic func;
	
		EquiRBFSO2() = default;
		
		explicit EquiRBFSO2( const VectorXd& centre ) : centre(centre) {}
		
		EquiRBFSO2( const VectorXd& weight, const VectorXd& centre ) : weight(weight), centre(centre) {}

		template<typename Derived>
		Matrix<typename Derived::Scalar,-1,1> operator()( const MatrixBase<Derived>& x ) const {
			typedef Derived::Scalar Scalar;
			Matrix<Scalar,2,1> Wtx = W.adjoint() * x;
			Matrix<double,2,1> Wtc = W.adjoint() * centre;
			Scalar rR = Wtx.squaredNorm() * Wtc.squaredNorm();
			if( rR == Scalar() ) {
				return 2.0 * M_PI * func( (x-centre.cast<Scalar>()).norm() ) * ( weight - W * ( W.adjoint() * weight ) ).cast<Scalar>();
			}
			Matrix<Scalar,2,2> rotation = ComputeRotation( Wtc, Wtx );
			return (-M_PI / sqrt(rR)) * ( W * ( rotation * ( W.adjoint() * weight ) ) );
		}

		MatrixXd Derivative( const VectorXd& x ) const {
			return AutoDerivative( *this, x );
		}

	protected:

		template<typename Scalar>
		static Matrix<Scalar,2,2> ComputeRotation( Scalar theta ) {
			Scalar sintheta = sin(theta);
			Scalar costheta = cos(theta);
			Matrix<Scalar,2,2> rotation;
			rotation(0,0) = costheta;
			rotation(1,0) = sintheta;
			rotation(0,1) = -sintheta;
			rotation(1,1) = costheta;
			return rotation;
		}

		template<typename Derived1,typename Derived2>
		static Matrix<typename internal::scalar_product_traits<typename Derived1::Scalar,typename Derived2::Scalar>::ReturnType,2,2> ComputeRotation( const MatrixBase<Derived1>& x, const MatrixBase<Derived2>& y ) {
			return ComputeRotation( acos( x.dot(y) / (x.norm()*y.norm()) ) );
		}
	};



}

#endif
