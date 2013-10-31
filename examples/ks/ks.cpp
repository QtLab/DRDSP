#include "ks.h"
#include <iostream>
#include <Eigen/LU>
#include <DRDSP/misc.h>

using namespace std;
using namespace Eigen;

void KSWrap::operator()( VectorXd& state ) const {
	for(uint32_t i=0;i<state.size()-1;i++)
		Wrap(state(i),-M_PI,M_PI);
}

// Embeddings

VectorXd FlatEmbedding::Evaluate( const VectorXd &state ) const {
	VectorXd X(eDim);
	
	uint32_t k=0;
	for(uint32_t i=0;i<oDim;i++) {
		X(k++) = cos(state(i));
		X(k++) = sin(state(i));
	}

	return std::move(X);
}

// Embedding First Derivatives

MatrixXd FlatEmbedding::Derivative( const VectorXd &state ) const {
	MatrixXd res;
	res.setZero(eDim,oDim);

	VectorXd X = Evaluate(state);

	uint32_t k=0;
	for(uint32_t i=0;i<oDim;i++) {
		res(k,i) = -X(k+1);
		res(k+1,i) = X(k);
		k += 2;
	}

	return std::move(res);
}

// Embedding First Derivatives Adjoint

MatrixXd FlatEmbedding::DerivativeAdjoint( const VectorXd &state ) const {
	return Derivative(state).transpose();
}

// Embedding Second Derivatives

MatrixXd FlatEmbedding::Derivative2( const VectorXd &x, uint32_t mu ) const {
	MatrixXd res;
	res.setZero(oDim,oDim);

	uint32_t k = mu / 2;

	if( IsOdd(mu) ) {
		res(k,k) = -sin(x(k));
	} else {
		res(k,k) = -cos(x(k));
	}

	return std::move(res);
}

// Model Functions

double KS::Forcing( const VectorXd& state ) const {
	return sin( state(dimension-1) );
}

double KS::ForcingDerivative( const VectorXd& state ) const {
	return cos( state(dimension-1) );
}

double KS::Phase( uint32_t i, const VectorXd& state ) const {
	return state(i);
}
	
complex<double> KS::ComplexOrderParameter( const VectorXd& state ) const {
	complex<double> z(0,0);
	for(uint32_t i=0;i<numOscillators;i++) {
		z += complex<double>( cos(Phase(i,state)), sin(Phase(i,state)) ) ;
	}
	z /= numOscillators;
	return std::move(z);
}

VectorXd KS::VectorField( const VectorXd& state, const VectorXd& parameter ) {
	complex<double> z = ComplexOrderParameter(state);
	VectorXd s(numOscillators);
	double psi = arg(z);
	for(uint32_t i=0;i<numOscillators;i++)
		s(i) = sin( psi - Phase(i,state) );
	VectorXd res(dimension);
	res.head(numOscillators) = frequencies + interactionStrengths * Forcing(state) + K * abs(z) * s;
	res(dimension-1) = forcingFrequency;
	return std::move(res);
}
	
MatrixXd KS::Partials( const VectorXd& state, const VectorXd& parameter ) {
	MatrixXd res;
	res.setZero(dimension,dimension);
			
	complex<double> z = ComplexOrderParameter(state);
	double psi = arg(z);
	double r = abs(z);
	double forcing = Forcing(state);
	for(uint32_t i=0;i<numOscillators;i++) {
		for(uint32_t j=0;j<numOscillators;j++) {
			res(i,j) = interactionStrengths(i) * forcing + K * MeanAmplitudeDerivative(j,state,psi) * sin( psi - Phase(i,state) ) + K  * cos( psi - Phase(i,state) ) * ( MeanPhaseDerivative(j,state,psi) - delta(i,j) * r );
		}
		res(i,dimension-1) = interactionStrengths(i) * ForcingDerivative(state);
	}
	return std::move(res);
}

double KS::delta( uint32_t i, uint32_t j ) const {
	return (i==j)?1.0:0.0;
}

double KS::MeanAmplitudeDerivative( uint32_t j, const VectorXd &state, double psi ) const {
	return ( cos(Phase(j,state)) * sin(psi) - sin(Phase(j,state)) * cos(psi) ) / numOscillators;
}

double KS::MeanPhaseDerivative( uint32_t j, const VectorXd &state, double psi ) const {
	return ( cos(Phase(j,state)) * cos(psi) + sin(Phase(j,state)) * sin(psi) ) / numOscillators;
}




