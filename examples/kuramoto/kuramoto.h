#ifndef INCLUDED_KURAMOTO
#define INCLUDED_KURAMOTO
#include <DRDSP/dynamics/model.h>
#include <complex>

using namespace std;
using namespace DRDSP;

struct KuramotoWrap : WrapFunction<VectorXd> {
	void operator()( VectorXd& x ) const;
};

struct Kuramoto : ModelParameterized {
	Kuramoto();
	explicit Kuramoto( uint32_t N );
	void Create( uint32_t N );
	VectorXd VectorField( const VectorXd& state, const VectorXd& parameter );
	MatrixXd Partials( const VectorXd& state, const VectorXd& parameter );

protected:
	KuramotoWrap wrap;
	VectorXd frequencies;
	VectorXd interactionStrengths;
	double K;
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

struct KuramotoFlat : ModelParameterizedEmbedded {
	explicit KuramotoFlat( uint32_t N ) : ModelParameterizedEmbedded(kuramoto,embedFlat), kuramoto(N), embedFlat(N+1) {}
protected:
	Kuramoto kuramoto;
	FlatEmbedding embedFlat;
};

#endif
