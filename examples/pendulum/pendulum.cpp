#include "pendulum.h"
#include <iostream>
#include <Eigen/LU>

using namespace std;
using namespace Eigen;

// Embeddings

VectorXd FlatEmbedding::Evaluate( const VectorXd &state ) const {
	VectorXd X(eDim);
	
	X(0) = cos(state(0));
	X(1) = sin(state(0));
	X(2) = cos(state(1));
	X(3) = sin(state(1));
	X(4) = cos(state(2));
	X(5) = sin(state(2));
	X(6) = state(3);
	X(7) = state(4);

	return std::move(X);
}

VectorXd DoughnutEmbedding::Evaluate( const VectorXd &state ) const {
	VectorXd y(eDim);

	double phi = state(0);
	double theta = state(1);
	double psi = state(2);
	double vp = state(3);
	double vt = state(4);

	double cosphi = cos(phi);
	double costheta = cos(theta);
	
	y(0) = (R2+(R1+cosphi)*costheta)*cos(psi);
	y(1) = (R2+(R1+cosphi)*costheta)*sin(psi);
	y(2) = (R1+cosphi)*sin(theta);
	y(3) = sin(phi);
	y(4) = vp;
	y(5) = vt;

	return std::move(y);
}

// Embedding First Derivatives

MatrixXd FlatEmbedding::Derivative( const VectorXd &state ) const {
	MatrixXd res;
	res.setZero(eDim,oDim);

	VectorXd X = Evaluate(state);

	res(0,0) = -X(1);
	res(1,0) = X(0);
	res(2,1) = -X(3);
	res(3,1) = X(2);
	res(4,2) = -X(5);
	res(5,2) = X(4);

	res(6,3) = 1;
	res(7,4) = 1;

	return std::move(res);
}

MatrixXd DoughnutEmbedding::Derivative( const VectorXd &x ) const {
	MatrixXd res;
	res.setZero(eDim,oDim);

	double cosphi = cos(x(0));
	double costheta = cos(x(1));
	double cospsi = cos(x(2));
	double sinphi = sin(x(0));
	double sintheta = sin(x(1));
	double sinpsi = sin(x(2));

	res(0,0) = -sinphi*costheta*cospsi;
	res(0,1) = -(R1+cosphi)*sintheta*cospsi;
	res(0,2) = -(R2+(R1+cosphi)*costheta)*sinpsi;
	
	res(1,0) = -sinphi*costheta*sinpsi;
	res(1,1) = -(R1+cosphi)*sintheta*sinpsi;
	res(1,2) = (R2+(R1+cosphi)*costheta)*cospsi;

	res(2,0) = -sinphi*sintheta;
	res(2,1) = (R1+cosphi)*costheta;

	res(3,0) = cosphi;

	res(4,3) = 1;
	res(5,4) = 1;
	
	return std::move(res);
}


// Embedding First Derivatives Adjoint

MatrixXd FlatEmbedding::DerivativeAdjoint( const VectorXd &state ) const {
	return Derivative(state).transpose();
}

MatrixXd DoughnutEmbedding::DerivativeAdjoint( const VectorXd &state ) const {
	return Derivative(state).transpose(); // check me
}

// Embedding Second Derivatives

MatrixXd FlatEmbedding::Derivative2( const VectorXd &x, uint32_t mu ) const {
	MatrixXd res;
	res.setZero(oDim,oDim);

	switch( mu ) {
		case 0:
			res(0,0) = -cos(x(0));
			break;
		case 1:
			res(0,0) = -sin(x(0));
			break;
		case 2:
			res(1,1) = -cos(x(1));
			break;
		case 3:
			res(1,1) = -sin(x(1));
			break;
		case 4:
			res(2,2) = -cos(x(2));
			break;
		case 5:
			res(2,2) = -sin(x(2));
			break;
	}

	return std::move(res);
}

MatrixXd DoughnutEmbedding::Derivative2( const VectorXd &x, uint32_t mu ) const {
	MatrixXd res;
	res.setZero(oDim,oDim);

	double cosphi = cos(x(0));
	double costheta = cos(x(1));
	double cospsi = cos(x(2));
	double sinphi = sin(x(0));
	double sintheta = sin(x(1));
	double sinpsi = sin(x(2));

	switch( mu ) {
		case 0:
			res(0,0) =				-cosphi*costheta*cospsi;
			res(1,0) = res(0,1) =	sinphi*sintheta*cospsi;
			res(2,0) = res(0,2) =	sinphi*costheta*sinpsi;
			res(1,1) =				-(R1+cosphi)*costheta*cospsi;
			res(2,1) = res(1,2) =	(R1+cosphi)*sintheta*sinpsi;
			res(2,2) =				-(R2+(R1+cosphi)*costheta)*cospsi;
			break;
		case 1:
			res(0,0) =				-cosphi*costheta*sinpsi;
			res(1,0) = res(0,1) =	sinphi*sintheta*sinpsi;
			res(2,0) = res(0,2) =	-sinphi*costheta*cospsi;
			res(1,1) =				-(R1+cosphi)*costheta*sinpsi;
			res(2,1) = res(1,2) =	-(R1+cosphi)*sintheta*cospsi;
			res(2,2) =				-(R2+(R1+cosphi)*costheta)*sinpsi;
			break;
		case 2:
			res(0,0) =				-cosphi*sintheta;
			res(1,0) = res(0,1) =	-sinphi*costheta;
			res(1,1) =				-(R1+cosphi)*sintheta;
			break;
		case 3:
			res(0,0) =				-sinphi;
			break;
	}

	return std::move(res);
}

// Model Functions

VectorXd Pendulum::VectorField( const VectorXd &state, const VectorXd &parameter ) {
	VectorXd res(dimension);
	res(0) = phiDot(state);
	res(1) = thetaDot(state);
	res(2) = psiDot(parameter);
	res(3) = vpDot(state,parameter);
	res(4) = vtDot(state,parameter);
	return std::move(res);
}

MatrixXd Pendulum::Partials( const VectorXd &state, const VectorXd &parameter ) {
	MatrixXd res;
	res.setZero(dimension,dimension);

	res(0,3) = 1;
	res(1,4) = 1;

	double sinphi = sin(state(0));
	double cosphi = cos(state(0));
	double sintheta = sin(state(1));
	double costheta = cos(state(1));

	double tf1 = f1(state(0));
	double tf1d = f1d(state(0));
	double tf2 = f2(state(0));
	double tf2d = f2d(state(0));
	double tf3 = f3(state(2),parameter);
	double tf3d = f3d(state(2),parameter);
	double tg1 = g1(state(0));
	double tg1d = g1d(state(0));
	double tf4 = f4(state(0));
	double tf4d = f4d(state(0));

	res(3,0) = -tg1d*state(4)*state(4)-tf3*cosphi*costheta;
	res(3,1) = tf3*sinphi*sintheta;
	res(3,2) = -tf3d*sinphi*costheta;
	res(3,3) = -delta2;
	res(3,4) = -2.0*tg1*state(4);

	res(4,0) = ( tf1d*state(3)*state(4) - tf2d*tf3*sintheta )/tf4 - (( tf1*state(3)*state(4)-delta1*state(4)-tf2*tf3*sintheta )*tf4d )/(tf4*tf4);
	res(4,1) = -tf2*tf3*costheta / tf4;
	res(4,2) = -tf2*tf3d*sintheta / tf4;
	res(4,3) = tf1*state(4) / tf4;
	res(4,4) = (tf1*state(3) - delta1) / tf4;

	return std::move(res);
}

double Pendulum::f1( double p ) const {
	double s = sin(p);
	double c = cos(p);
	return (3.0*length+2.0*c)*s;
}

double Pendulum::f2( double p ) const {
	return (mass+2.0)*length+cos(p);
}

double Pendulum::f3( double p, const VectorXd &b ) const {
	double omega = b(0);
	return 1.0 + A*omega*omega*cos(p);
}

double Pendulum::f4( double p ) const {
	double c = cos(p);
	return (mass+3.0)*length*length + 3.0*length*c + c*c;
}

inline double Pendulum::g1( double p ) const {
	return f1(p)*0.5;
}

inline double Pendulum::phiDot( const VectorXd &x ) const {
	return x(3);
}

inline double Pendulum::thetaDot( const VectorXd &x ) const {
	return x(4);
}


inline double Pendulum::psiDot( const VectorXd &b ) const {
	return b(0);
}

double Pendulum::vpDot( const VectorXd &x, const VectorXd &b ) const {
	return -g1(x[0])*x[4]*x[4] - delta2*x[3] - f3(x[2],b)*sin(x[0])*cos(x[1]);
}

double Pendulum::vtDot( const VectorXd &x, const VectorXd &b) const {
	return (f1(x[0])*x[4]*x[3] - delta1*x[4] - f2(x[0])*f3(x[2],b)*sin(x[1]) )/f4(x[0]);
}

double Pendulum::f1d( double p ) const {
	double s = sin(p);
	double c = cos(p);
	return (3.0*length+2.0*c)*c - 2.0*s*s;
}

inline double Pendulum::f2d( double p ) const {
	return -sin(p);
}

double Pendulum::f3d( double p, const VectorXd &b ) const {
	double omega = b(0);
	//double A = b(1);
	return -A*omega*omega*sin(p);
}

inline double Pendulum::f4d( double p ) const {
	return -f1(p);
}

inline double Pendulum::g1d( double p ) const {
	return f1d(p)*0.5;
}




