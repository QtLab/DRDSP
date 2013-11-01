#ifndef INCLUDED_KS
#define INCLUDED_KS
#include <DRDSP/dynamics/model.h>
#include <complex>

using namespace std;
using namespace DRDSP;

struct KSWrap : WrapFunction<VectorXd> {
	void operator()( VectorXd& x ) const;
};

struct KS : ModelParameterized {
	KS();
	explicit KS( uint32_t N );
	void Create( uint32_t N );
	VectorXd VectorField( const VectorXd& state, const VectorXd& parameter );
	MatrixXd Partials( const VectorXd& state, const VectorXd& parameter );

protected:
	KSWrap wrap;
	VectorXd frequencies;
	VectorXd interactionStrengths;
	double forcingFrequency;
	uint32_t numOscillators;

	double Forcing( const VectorXd& state ) const;
	double ForcingDerivative( const VectorXd& state ) const;
	double Phase( uint32_t i, const VectorXd& state ) const;
	complex<double> ComplexOrderParameter( const VectorXd& state ) const;
	double MeanAmplitudeDerivative( uint32_t j, const VectorXd &state, double psi ) const;
	double MeanPhaseDerivative( uint32_t j, const VectorXd &state, double psi ) const;
};

struct FlatEmbedding : Embedding {
	explicit FlatEmbedding( uint32_t n ) : Embedding(n,2*n) {}
	VectorXd Evaluate( const VectorXd &x ) const;
	MatrixXd Derivative( const VectorXd &x ) const;
	MatrixXd DerivativeAdjoint( const VectorXd &x ) const;
	MatrixXd Derivative2( const VectorXd &x, uint32_t mu ) const;
};

struct KSFlat : ModelParameterizedEmbedded {
	explicit KSFlat( uint32_t N ) : ModelParameterizedEmbedded(ks,embedFlat), ks(N), embedFlat(N+1) {}
protected:
	KS ks;
	FlatEmbedding embedFlat;
};

#endif
