#ifndef INCLUDED_PENDULUM
#define INCLUDED_PENDULUM
#include <DRDSP/types.h>
#include <DRDSP/dynamics/model_orig.h>

using namespace DRDSP;

struct Pendulum : ModelParameterized {

	Pendulum() : ModelParameterized(5,1), length(1.1), mass(7.0), A(0.15), delta1(0.245), delta2(0.245) {}
	VectorXd VectorField( const VectorXd &state, const VectorXd &parameter );
	MatrixXd Partials( const VectorXd &state, const VectorXd &parameter );

protected:
	// fixed parameters
	double length, mass, A, delta1, delta2;

	double f1( double p );
	double f2( double p );
	double f3( double p, const VectorXd &b );
	double f4( double p );
	double g1( double p );

	double phiDot( const VectorXd &x );
	double thetaDot( const VectorXd &x );
	double psiDot( const VectorXd &b );
	double vpDot( const VectorXd &x, const VectorXd &b );
	double vtDot( const VectorXd &x, const VectorXd &b);

	VectorXd G( const VectorXd &theta, const VectorXd &b );
	MatrixXd DG( const VectorXd &x, const VectorXd &b );

	double f1d( double p );
	double f2d( double p );
	double f3d( double p, const VectorXd &b );
	double f4d( double p );
	double g1d( double p );

};

struct FlatEmbedding : Embedding {
	FlatEmbedding() : Embedding(5,8) {}
	VectorXd Evaluate( const VectorXd &x ) const;
	MatrixXd Derivative( const VectorXd &x ) const;
	MatrixXd DerivativeAdjoint( const VectorXd &x ) const;
	MatrixXd Derivative2( const VectorXd &x, uint32_t mu ) const;
	
};

struct DoughnutEmbedding : Embedding {
	double R1, R2;
	DoughnutEmbedding() : Embedding(5,6), R1(2.0), R2(4.0) {}
	VectorXd Evaluate( const VectorXd &x ) const;
	MatrixXd Derivative( const VectorXd &x ) const;
	MatrixXd DerivativeAdjoint( const VectorXd &x ) const;
	MatrixXd Derivative2( const VectorXd &x, uint32_t mu ) const;
};

struct PendulumFlat : ModelParameterizedEmbedded {

	PendulumFlat() : ModelParameterizedEmbedded(pendulum,embedFlat) {}

protected:
	Pendulum pendulum;
	FlatEmbedding embedFlat;
};

struct PendulumDoughnut : ModelParameterizedEmbedded {

	PendulumDoughnut() : ModelParameterizedEmbedded(pendulum,embedDoughnut) {}

protected:
	Pendulum pendulum;
	DoughnutEmbedding embedDoughnut;
};

#endif
