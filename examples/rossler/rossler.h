#ifndef INCLUDED_ROSSLER
#define INCLUDED_ROSSLER
#include <DRDSP/dynamics/model.h>

using namespace DRDSP;

struct RosslerJacobian {
	double a,b,c;
	RosslerJacobian() : a(0.1), b(0.1), c(4.0) {}
	explicit RosslerJacobian( double c ) : a(0.1), b(0.1), c(c) {}
	Matrix3d operator()( const Vector3d& x ) const;
};

struct RosslerModel : Model<Vector3d> {
	typedef double Time;

	double a,b,c;

	RosslerModel() : Model<Vector3d>(3), a(0.1), b(0.1), c(4.0) {}
	explicit RosslerModel( double c ) : Model<Vector3d>(3), a(0.1), b(0.1), c(c) {}
	Vector3d operator()( const Vector3d& x ) const;
	Matrix3d Partials( const Vector3d& x ) const;
	RosslerJacobian Jacobian() const { return RosslerJacobian(c); }
};

struct RosslerFamily : Family<RosslerModel> {
	
	RosslerFamily() : Family<RosslerModel>(3,1) {}
	
	RosslerModel operator()( const VectorXd& parameter ) const {
		return RosslerModel( parameter[0] );
	}

};

#endif
