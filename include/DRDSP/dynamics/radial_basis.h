#ifndef INCLUDED_DYNAMICS_RADIAL_BASIS
#define INCLUDED_DYNAMICS_RADIAL_BASIS
#include "../types.h"

namespace DRDSP {

	struct ThinPlateSpline {
		double operator()( double r ) const;
		double Derivative( double r ) const;
	};

	struct PolyharmonicSpline3 {
		double operator()( double r ) const;
		double Derivative( double r ) const;
	};

	struct Gaussian  {
		double scale;
		Gaussian();
		explicit Gaussian( double scale );
		double operator()( double r ) const;
		double Derivative( double r ) const;
	};

	struct Multiquadratic {
		double scale;
		Multiquadratic();
		explicit Multiquadratic( double scale );
		double operator()( double r ) const;
		double Derivative( double r ) const;
	};

	struct InverseQuadratic {
		double scale;
		InverseQuadratic();
		explicit InverseQuadratic( double scale );
		double operator()( double r ) const;
		double Derivative( double r ) const;
	};

	struct InverseMultiquadratic {
		double scale;
		InverseMultiquadratic();
		explicit InverseMultiquadratic( double scale );
		double operator()( double r ) const;
		double Derivative( double r ) const;
	};

	template<typename F>
	struct RadialFunction {
		VectorXd centre;
		F f;
	
		RadialFunction() = default;
		
		explicit RadialFunction( const F& f ) : f(f) {}
		
		explicit RadialFunction( const VectorXd& centre ) : centre(centre) {}
		
		RadialFunction( const F& f, const VectorXd& centre ) : f(f), centre(centre) {}

		double operator()( const VectorXd& x ) const {
			return f( (x-centre).norm() );
		}

		VectorXd Derivative( const VectorXd& x ) const {
			VectorXd r = x - centre;
			double rnorm = r.norm();
			if( rnorm == 0.0 ) return VectorXd::Zero(x.size());
			return ( f.Derivative( rnorm ) / rnorm ) * r;
		}

	};

}

#endif
