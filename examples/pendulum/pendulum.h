#ifndef INCLUDED_PENDULUM
#define INCLUDED_PENDULUM
#include <DRDSP/types.h>
#include <DRDSP/dynamics/model.h>
#include <DRDSP/dynamics/dynamicalSystem.h>

using namespace DRDSP;

struct PendulumWrap {
	void operator()( VectorXd& x ) const;
};

struct Pendulum : Model<> {
	PendulumWrap pendulumWrap;

	Pendulum() : Pendulum(1.8) {}

	explicit Pendulum( double Omega ) : Model<>(5), Omega(Omega), length(1.1), mass(7.0), A(0.15), delta1(0.245), delta2(0.245) {}

	VectorXd operator()( const VectorXd& state ) const;
	MatrixXd Partials( const VectorXd& state ) const;

protected:
	// parameters
	double Omega, length, mass, A, delta1, delta2;

	double f1( double p ) const;
	double f2( double p ) const;
	double f3( double p ) const;
	double f4( double p ) const;
	double g1( double p ) const;

	double phiDot( const VectorXd& x ) const;
	double thetaDot( const VectorXd& x ) const;
	double psiDot() const;
	double vpDot( const VectorXd& x ) const;
	double vtDot( const VectorXd& x ) const;

	VectorXd G( const VectorXd& theta ) const;
	MatrixXd DG( const VectorXd& x ) const;

	double f1d( double p ) const;
	double f2d( double p ) const;
	double f3d( double p ) const;
	double f4d( double p ) const;
	double g1d( double p ) const;

};

struct PendulumFamily : Family<Pendulum> {

	PendulumFamily() : Family<Pendulum>(5,1) {}

	Pendulum operator()( const VectorXd& parameter ) const {
		return Pendulum( parameter[0] );
	}
};

struct FlatEmbedding : Embedding {
	FlatEmbedding() : Embedding(5,8) {}
	VectorXd operator()( const VectorXd& x ) const;
	MatrixXd Derivative( const VectorXd& x ) const;
	MatrixXd DerivativeAdjoint( const VectorXd& x ) const;
	MatrixXd Derivative2( const VectorXd& x, uint32_t mu ) const;
	
};

struct DoughnutEmbedding : Embedding {
	double R1, R2;
	DoughnutEmbedding() : Embedding(5,6), R1(2.0), R2(4.0) {}
	VectorXd operator()( const VectorXd& x ) const;
	MatrixXd Derivative( const VectorXd& x ) const;
	MatrixXd DerivativeAdjoint( const VectorXd& x ) const;
	MatrixXd Derivative2( const VectorXd& x, uint32_t mu ) const;
};

typedef RKDynamicalSystem<SolverFunctionFromModel<Pendulum>,PendulumWrap> PendulumSolver;

#endif
