#ifndef INCLUDED_ROSSLER
#define INCLUDED_ROSSLER
#include <DRDSP/dynamics/model.h>
#include <DRDSP/auto_diff.h>

using namespace DRDSP;

struct RosslerModel : Model<Vector3d> {
	double a,b,c;

	RosslerModel() : Model<Vector3d>(3), a(0.1), b(0.1), c(4.0) {}
	
	explicit RosslerModel( double c ) : Model<Vector3d>(3), a(0.1), b(0.1), c(c) {}

	template<typename Derived>
	Matrix<typename Derived::Scalar,3,1> operator()( const MatrixBase<Derived>& x ) const {
		Matrix<typename Derived::Scalar,3,1> v;
		v[0] = -x[1] - x[2];
		v[1] = x[0] + a * x[1];
		v[2] = b + x[2] * ( x[0] - c );
		return v;
	}

	Matrix3d Partials( const Vector3d& x ) const {
		return AutoDerivative( *this, x );
	}
};

struct RosslerFamily : Family<RosslerModel> {
	
	RosslerFamily() : Family<RosslerModel>(3,1) {}
	
	RosslerModel operator()( const VectorXd& parameter ) const {
		return RosslerModel( parameter[0] );
	}

};

#endif
