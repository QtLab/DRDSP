#pragma warning( disable : 4530 )
#include "rossler.h"

Vector3d RosslerModel::operator()( const Vector3d& x ) const {
	VectorXd v(3);
	v(0) = -x(1) - x(2);
	v(1) = x(0) + a * x(1);
	v(2) = b + x(2) * ( x(0) - c );
	return v;
}

Matrix3d RosslerModel::Partials( const Vector3d& x ) const {
	Matrix3d res;

	res(0,0) = 0.0;
	res(0,1) = -1.0;
	res(0,2) = -1.0;
	res(1,0) = 1.0;
	res(1,1) = a;
	res(1,2) = 0.0;
	res(2,0) = x(2);
	res(2,1) = 0.0;
	res(2,2) = x(0) - c;

	return res;
}

Matrix3d RosslerJacobian::operator()( const Vector3d& x ) const {
	Matrix3d res;

	res(0,0) = 0.0;
	res(0,1) = -1.0;
	res(0,2) = -1.0;
	res(1,0) = 1.0;
	res(1,1) = a;
	res(1,2) = 0.0;
	res(2,0) = x(2);
	res(2,1) = 0.0;
	res(2,2) = x(0) - c;

	return res;
}
