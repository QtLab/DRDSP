#ifndef INCLUDED_DYNAMICS_RADIAL_BASIS
#define INCLUDED_DYNAMICS_RADIAL_BASIS
#include "../types.h"

namespace DRDSP {

	struct Function {
		virtual double operator()( double r ) const = 0;
		virtual double Derivative( double r ) const = 0;
	};

	struct ThinPlateSpline : Function {
		double operator()( double r ) const;
		double Derivative( double r ) const;
	};

	struct PolyharmonicSpline3 : Function {
		double operator()( double r ) const;
		double Derivative( double r ) const;
	};

	struct Gaussian : Function {
		double scale;
		Gaussian();
		double operator()( double r ) const;
		double Derivative( double r ) const;
	};

	struct Multiquadratic : Function {
		double scale;
		Multiquadratic();
		double operator()( double r ) const;
		double Derivative( double r ) const;
	};

	struct InverseQuadratic : Function {
		double scale;
		InverseQuadratic();
		double operator()( double r ) const;
		double Derivative( double r ) const;
	};

	struct InverseMultiquadratic : Function {
		double scale;
		InverseMultiquadratic();
		double operator()( double r ) const;
		double Derivative( double r ) const;
	};

	struct RadialFunction {
		const Function* function;
		VectorXd centre;
	
		RadialFunction();
		RadialFunction( const Function& f );
		double operator()( const VectorXd& x ) const;
		VectorXd Derivative( const VectorXd& x ) const;
	protected:
		static ThinPlateSpline thinPlateSpline;
	};

}

#endif
