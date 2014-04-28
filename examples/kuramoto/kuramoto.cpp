#include "kuramoto.h"
#include <DRDSP/misc.h>
#include <cmath>

using namespace std;

VectorXd FlatEmbedding::operator()( const VectorXd &state ) const {
	VectorXd X(eDim);
	
	uint32_t k=0;
	for(uint32_t i=0;i<oDim;i++) {
		X(k++) = cos(state(i));
		X(k++) = sin(state(i));
	}

	return X;
}

MatrixXd FlatEmbedding::Derivative( const VectorXd &state ) const {
	MatrixXd res;
	res.setZero(eDim,oDim);

	VectorXd X = (*this)(state);

	uint32_t k=0;
	for(uint32_t i=0;i<oDim;i++) {
		res(k,i) = -X(k+1);
		res(k+1,i) = X(k);
		k += 2;
	}

	return res;
}

MatrixXd FlatEmbedding::DerivativeAdjoint( const VectorXd &state ) const {
	return Derivative(state).transpose();
}

MatrixXd FlatEmbedding::Derivative2( const VectorXd &x, uint32_t mu ) const {
	MatrixXd res;
	res.setZero(oDim,oDim);

	uint32_t k = mu / 2;

	if( IsOdd(mu) ) {
		res(k,k) = -sin(x(k));
	} else {
		res(k,k) = -cos(x(k));
	}

	return res;
}

void KuramotoWrap::operator()( VectorXd& state ) const {
	for(int i=0;i<state.size();i++)
		Wrap(state(i),-M_PI,M_PI);
}



double KuramotoBase::Forcing( const VectorXd& state ) const {
	return sin( state(numOscillators) );
}

double KuramotoBase::ForcingDerivative( const VectorXd& state ) const {
	return cos( state(numOscillators) );
}

double KuramotoBase::Phase( uint32_t i, const VectorXd& state ) const {
	return state(i);
}
	
complex<double> KuramotoBase::ComplexOrderParameter( const VectorXd& state ) const {
	complex<double> z(0,0);
	for(uint32_t i=0;i<numOscillators;i++) {
		z += complex<double>( cos(Phase(i,state)), sin(Phase(i,state)) ) ;
	}
	z /= numOscillators;
	return z;
}

double KuramotoBase::MeanAmplitudeDerivative( uint32_t j, const VectorXd &state, double psi ) const {
	return ( cos(Phase(j,state)) * sin(psi) - sin(Phase(j,state)) * cos(psi) ) / numOscillators;
}

double KuramotoBase::MeanPhaseDerivative( uint32_t j, const VectorXd &state, double psi ) const {
	return ( cos(Phase(j,state)) * cos(psi) + sin(Phase(j,state)) * sin(psi) ) / numOscillators;
}


VectorXd KuramotoA::operator()( const VectorXd& state ) const {
	complex<double> z = ComplexOrderParameter(state);
	VectorXd s(numOscillators);
	double psi = arg(z);
	for(uint32_t i=0;i<numOscillators;i++)
		s(i) = sin( psi - Phase(i,state) );
	VectorXd res(dimension);
	res.head(numOscillators) = frequencies + interactionStrengths * Forcing(state) + K * abs(z) * s;
	res(dimension-1) = forcingFrequency;
	return res;
}

MatrixXd KuramotoA::Partials( const VectorXd& state ) const {
	MatrixXd res;
	res.setZero(dimension,dimension);
			
	complex<double> z = ComplexOrderParameter(state);
	double psi = arg(z);
	double r = abs(z);
	double forcing = Forcing(state);
	for(uint32_t i=0;i<numOscillators;i++) {
		for(uint32_t j=0;j<numOscillators;j++) {
			res(i,j) = interactionStrengths(i) * forcing + K * MeanAmplitudeDerivative(j,state,psi) * sin( psi - Phase(i,state) ) + K  * cos( psi - Phase(i,state) ) * ( MeanPhaseDerivative(j,state,psi) - Delta<double>(i,j) * r );
		}
		res(i,dimension-1) = interactionStrengths(i) * ForcingDerivative(state);
	}
	return res;
}

VectorXd KuramotoB::operator()( const VectorXd& state ) const {
	complex<double> z = ComplexOrderParameter(state);
	VectorXd s(numOscillators);
	double psi = arg(z);
	double r = abs(z);
	for(uint32_t i=0;i<numOscillators;i++)
		s(i) = r * sin( psi - Phase(i,state) );
	VectorXd res(dimension);
	res.head(numOscillators) = frequencies + (interactionStrengths * Forcing(state) + K * VectorXd::Ones(numOscillators) ).cwiseProduct(s);
	res(dimension-1) = forcingFrequency;
	return res;
}

MatrixXd KuramotoB::Partials( const VectorXd& state ) const {
	MatrixXd res;
	res.setZero(dimension,dimension);
			
	complex<double> z = ComplexOrderParameter(state);
	double psi = arg(z);
	double r = abs(z);
	double forcing = Forcing(state);
	for(uint32_t i=0;i<numOscillators;i++) {
		double sinPsiThetai = sin( psi - Phase(i,state) );
		double cosPsiThetai = cos( psi - Phase(i,state) );
		for(uint32_t j=0;j<numOscillators;j++) {
			res(i,j) = (interactionStrengths(i) * forcing + K ) * ( MeanAmplitudeDerivative(j,state,psi) * sinPsiThetai + cosPsiThetai * ( MeanPhaseDerivative(j,state,psi) - Delta<double>(i,j) * r ) );
		}
		res(i,dimension-1) = interactionStrengths(i) * ForcingDerivative(state) * r * sinPsiThetai;
	}
	return res;
}
