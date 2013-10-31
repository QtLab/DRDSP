#ifndef INCLUDED_KS
#define INCLUDED_KS
#include <DRDSP/types.h>
#include <DRDSP/dynamics/model.h>
#include <complex>

using namespace std;
using namespace DRDSP;

struct KSWrap : WrapFunction<VectorXd> {
	void operator()( VectorXd& x ) const;
};

struct KS : ModelParameterized {
	KSWrap wrap;
	VectorXd frequencies;
	VectorXd interactionStrengths;
	double K, forcingFrequency;
	uint32_t numOscillators;

	double Forcing( const VectorXd& state ) const;
	double ForcingDerivative( const VectorXd& state ) const;
	double Phase( uint32_t i, const VectorXd& state ) const;
	complex<double> ComplexOrderParameter( const VectorXd& state ) const;
	VectorXd VectorField( const VectorXd& state, const VectorXd& parameter );
	MatrixXd Partials( const VectorXd& state, const VectorXd& parameter );
	double delta( uint32_t i, uint32_t j ) const;
	double MeanAmplitudeDerivative( uint32_t j, const VectorXd &state, double psi ) const;
	double MeanPhaseDerivative( uint32_t j, const VectorXd &state, double psi ) const;

	KS() : ModelParameterized(&wrap), numOscillators(0), K(1.0), forcingFrequency(1.0) {}
	KS( uint32_t N ) : ModelParameterized(&wrap,N+1,1), numOscillators(N), frequencies(N), interactionStrengths(N), K(1.0), forcingFrequency(1.0) {}

	void Create( uint32_t N ) {
		frequencies.setZero(N);
		interactionStrengths.setZero(N);
		numOscillators = N;
	}

protected:

};

struct FlatEmbedding : Embedding {
	FlatEmbedding() : Embedding(5,8) {}
	VectorXd Evaluate( const VectorXd &x ) const;
	MatrixXd Derivative( const VectorXd &x ) const;
	MatrixXd DerivativeAdjoint( const VectorXd &x ) const;
	MatrixXd Derivative2( const VectorXd &x, uint32_t mu ) const;
	
};


struct KSFlat : ModelParameterizedEmbedded {

	KSFlat() : ModelParameterizedEmbedded(ks,embedFlat) {}

protected:
	KS ks;
	FlatEmbedding embedFlat;
};


#endif
