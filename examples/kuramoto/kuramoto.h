#ifndef INCLUDED_KURAMOTO
#define INCLUDED_KURAMOTO
#include <DRDSP/dynamics/model.h>
#include <complex>

using namespace std;
using namespace DRDSP;

struct KuramotoWrap {
	void operator()( VectorXd& x ) const;
};

struct KuramotoBase : Model<> {
	KuramotoWrap wrap;
	VectorXd frequencies;
	VectorXd interactionStrengths;
	double K, forcingFrequency;
	uint32_t numOscillators;

	KuramotoBase() : KuramotoBase(1) {}

	explicit KuramotoBase( uint32_t N ) :
		Model<>(N+1),
		numOscillators(N),
		frequencies(N),
		interactionStrengths(N),
		K(1.0),
		forcingFrequency(1.0)
	{}

protected:
	double Forcing( const VectorXd& state ) const;
	double ForcingDerivative( const VectorXd& state ) const;
	double Phase( uint32_t i, const VectorXd& state ) const;
	complex<double> ComplexOrderParameter( const VectorXd& state ) const;
	double MeanAmplitudeDerivative( uint32_t j, const VectorXd& state, double psi ) const;
	double MeanPhaseDerivative( uint32_t j, const VectorXd& state, double psi ) const;
};

struct KuramotoA : KuramotoBase {
	KuramotoA() = default;
	explicit KuramotoA( uint32_t N ) : KuramotoBase(N) {}
	VectorXd operator()( const VectorXd& state ) const;
	MatrixXd Partials( const VectorXd& state ) const;
};

struct KuramotoB : KuramotoBase {
	KuramotoB() = default;
	explicit KuramotoB( uint32_t N ) : KuramotoBase(N) {}
	VectorXd operator()( const VectorXd& state ) const;
	MatrixXd Partials( const VectorXd& state ) const;
};

struct FlatEmbedding : Embedding {
	explicit FlatEmbedding( uint32_t n ) : Embedding(n,2*n) {}
	VectorXd operator()( const VectorXd& x ) const;
	MatrixXd Derivative( const VectorXd& x ) const;
	MatrixXd DerivativeAdjoint( const VectorXd& x ) const;
	MatrixXd Derivative2( const VectorXd& x, uint32_t mu ) const;
};

struct KuramotoAFamily : Family<KuramotoA> {
	uint32_t N;
	
	explicit KuramotoAFamily( uint32_t N ) : Family<KuramotoA>(N+1,1), N(N) {}
	
	KuramotoA operator()( const VectorXd& parameter ) const {
		KuramotoA model(N);
		model.forcingFrequency = parameter[0];
		return model;
	}
};

struct KuramotoBFamily : Family<KuramotoB> {
	uint32_t N;
	
	explicit KuramotoBFamily( uint32_t N ) : Family<KuramotoB>(N+1,1), N(N) {}
	
	KuramotoB operator()( const VectorXd& parameter ) const {
		KuramotoB model(N);
		model.forcingFrequency = parameter[0];
		return model;
	}
};

typedef RKDynamicalSystem<KuramotoA,KuramotoWrap> KuramotoASolver;
typedef RKDynamicalSystem<KuramotoB,KuramotoWrap> KuramotoBSolver;

#endif
