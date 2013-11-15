#include <DRDSP/dynamics/radial_basis.h>
#include <cmath>

using namespace DRDSP;

ThinPlateSpline RadialFunction::thinPlateSpline;

RadialFunction::RadialFunction() : function(&thinPlateSpline) {
}

RadialFunction::RadialFunction( const Function& f ) : function(&f) {
}
	
double RadialFunction::operator()( const VectorXd& x ) const {
	return (*function)( (x-centre).norm() );
}

VectorXd RadialFunction::Derivative( const VectorXd& x ) const {
	VectorXd r = x - centre;
	double rnorm = r.norm();
	if( rnorm == 0.0 ) VectorXd::Zero(x.size());
	return ( function->Derivative( rnorm ) / rnorm ) * r;
}

double ThinPlateSpline::operator()( double r ) const {
	return r*r*log(r);
}

double ThinPlateSpline::Derivative( double r ) const {
	return r * ( 1.0 + 2.0 * log(r) );
}

double PolyharmonicSpline3::operator()( double r ) const {
	return r*r*r;
}

double PolyharmonicSpline3::Derivative( double r ) const {
	return 3.0*r*r;
}

Gaussian::Gaussian() : scale(1.0) {
}

double Gaussian::operator()( double r ) const {
	double x = scale * r;
	return exp(-x*x);
}

double Gaussian::Derivative( double r ) const {
	double x = scale * r;
	return -2.0 * scale * x * exp(-x*x);
}

Multiquadratic::Multiquadratic() : scale(1.0) {
}

double Multiquadratic::operator()( double r ) const {
	double x = scale * r;
	return sqrt( 1.0 + x*x );
}

double Multiquadratic::Derivative( double r ) const {
	double x = scale * r;
	return ( scale * x ) / sqrt( 1.0 + x*x );
}

InverseQuadratic::InverseQuadratic() : scale(1.0) {
}

double InverseQuadratic::operator()( double r ) const {
	double x = scale * r;
	return 1.0 / ( 1.0 + x*x );
}

double InverseQuadratic::Derivative( double r ) const {
	double x = scale * r;
	double y = 1.0 + x*x;
	return ( -2.0 * scale * x ) / ( y*y );
}

InverseMultiquadratic::InverseMultiquadratic() : scale(1.0) {
}

double InverseMultiquadratic::operator()( double r ) const {
	double x = scale * r;
	return 1.0 / sqrt( 1.0 + x*x );
}

double InverseMultiquadratic::Derivative( double r ) const {
	double x = scale * r;
	double y = 1.0 + x*x;
	return ( -scale * x ) / sqrt( y*y*y );
}



