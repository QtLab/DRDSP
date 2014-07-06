#include "pendulum.h"
#include <DRDSP/misc.h>

using namespace std;

void PendulumWrap::operator()( VectorXd& state ) const {
	Wrap(state(0),-M_PI,M_PI);
	Wrap(state(1),-M_PI,M_PI);
	Wrap(state(2),-M_PI,M_PI);
}

// Embeddings

VectorXd FlatEmbedding::operator()( const VectorXd& state ) const {
	VectorXd X(embedDim);
	
	X(0) = cos(state(0));
	X(1) = sin(state(0));
	X(2) = cos(state(1));
	X(3) = sin(state(1));
	X(4) = cos(state(2));
	X(5) = sin(state(2));
	X(6) = state(3);
	X(7) = state(4);

	return X;
}

VectorXd DoughnutEmbedding::operator()( const VectorXd& state ) const {
	VectorXd y(embedDim);

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

	return y;
}

// Embedding First Derivatives

MatrixXd FlatEmbedding::Derivative( const VectorXd& state ) const {
	MatrixXd res;
	res.setZero(embedDim,sourceDim);

	VectorXd X = (*this)(state);

	res(0,0) = -X(1);
	res(1,0) = X(0);
	res(2,1) = -X(3);
	res(3,1) = X(2);
	res(4,2) = -X(5);
	res(5,2) = X(4);

	res(6,3) = 1;
	res(7,4) = 1;

	return res;
}

MatrixXd DoughnutEmbedding::Derivative( const VectorXd& x ) const {
	MatrixXd res;
	res.setZero(embedDim,sourceDim);

	double sinphi = sin(x(0));
	double cosphi = cos(x(0));
	double sintheta = sin(x(1));
	double costheta = cos(x(1));
	double sinpsi = sin(x(2));
	double cospsi = cos(x(2));

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
	
	return res;
}


// Embedding First Derivatives Adjoint

MatrixXd FlatEmbedding::DerivativeAdjoint( const VectorXd& state ) const {
	return Derivative(state).transpose();
}

MatrixXd DoughnutEmbedding::DerivativeAdjoint( const VectorXd& state ) const {
	return Derivative(state).transpose(); // check me
}

// Embedding Second Derivatives

MatrixXd FlatEmbedding::Derivative2( const VectorXd& x, uint32_t mu ) const {
	MatrixXd res;
	res.setZero(sourceDim,sourceDim);

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

	return res;
}

MatrixXd DoughnutEmbedding::Derivative2( const VectorXd& x, uint32_t mu ) const {
	MatrixXd res;
	res.setZero(sourceDim,sourceDim);

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

	return res;
}

// Model Functions

VectorXd Pendulum::operator()( const VectorXd& state ) const {
	double sinphi = sin(state(0));
	double cosphi = cos(state(0));
	double sintheta = sin(state(1));
	double costheta = cos(state(1));
	double cospsi = cos(state(2));
	
	VectorXd res(stateDim);
	res(0) = state[3];
	res(1) = state[4];
	res(2) = Omega;
	res(3) = vpDot( sinphi, cosphi, costheta, cospsi, state[3], state[4] );
	res(4) = vtDot( sinphi, cosphi, sintheta, cospsi, state[3], state[4] );
	return res;
}

MatrixXd Pendulum::Partials( const VectorXd& state ) const {
	MatrixXd res;
	res.setZero(stateDim,stateDim);

	res(0,3) = 1;
	res(1,4) = 1;

	double sinphi = sin(state(0));
	double cosphi = cos(state(0));
	double sintheta = sin(state(1));
	double costheta = cos(state(1));
	double sinpsi = sin(state(2));
	double cospsi = cos(state(2));

	double tf1 = f1(sinphi,cosphi);
	double tf1d = f1d(sinphi,cosphi);
	double tf2 = f2(cosphi);
	double tf2d = -sinphi;
	double tf3 = f3(cospsi);
	double tf3d = f3d(sinpsi);
	double tf4 = f4(cosphi);
	double tf4d = -tf1;

	res(3,0) = -0.5*tf1d*state(4)*state(4)-tf3*cosphi*costheta;
	res(3,1) = tf3*sinphi*sintheta;
	res(3,2) = -tf3d*sinphi*costheta;
	res(3,3) = -delta2;
	res(3,4) = -tf1*state(4);

	res(4,0) = ( tf1d*state(3)*state(4) - tf2d*tf3*sintheta )/tf4 - (( tf1*state(3)*state(4)-delta1*state(4)-tf2*tf3*sintheta )*tf4d )/(tf4*tf4);
	res(4,1) = -tf2*tf3*costheta / tf4;
	res(4,2) = -tf2*tf3d*sintheta / tf4;
	res(4,3) = tf1*state(4) / tf4;
	res(4,4) = (tf1*state(3) - delta1) / tf4;

	return res;
}

double Pendulum::f1( double sinphi, double cosphi ) const {
	return (3.0*length+2.0*cosphi)*sinphi;
}

double Pendulum::f2( double cosphi ) const {
	return (mass+2.0)*length+cosphi;
}

double Pendulum::f3( double cospsi ) const {
	return 1.0 + A*Omega*Omega*cospsi;
}

double Pendulum::f4( double cosphi ) const {
	return (mass+3.0)*length*length + 3.0*length*cosphi + cosphi*cosphi;
}

double Pendulum::vpDot( double sinphi, double cosphi, double costheta, double cospsi, double vp, double vt ) const {
	return -0.5*f1(sinphi,cosphi)*vt*vt - delta2*vp - f3(cospsi)*sinphi*costheta;
}

double Pendulum::vtDot( double sinphi, double cosphi, double sintheta, double cospsi, double vp, double vt ) const {
	return (f1(sinphi,cosphi)*vt*vp - delta1*vt - f2(cosphi)*f3(cospsi)*sintheta )/f4(cosphi);
}

double Pendulum::f1d( double sinphi, double cosphi ) const {
	return (3.0*length+2.0*cosphi)*cosphi - 2.0*sinphi*sinphi;
}

double Pendulum::f3d( double sinpsi ) const {
	return -A*Omega*Omega*sinpsi;
}


