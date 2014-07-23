#ifndef INCLUDED_PENDULUM
#define INCLUDED_PENDULUM
#include <DRDSP/types.h>
#include <DRDSP/dynamics/model.h>

using namespace DRDSP;

struct PendulumWrap {
	void operator()( VectorXd& x ) const;
};

struct Pendulum : Model<> {

	Pendulum() : Model<>(5) {}

	VectorXd operator()( const VectorXd& state ) const;
	MatrixXd Partials( const VectorXd& state ) const;

	friend struct PendulumFamily;

protected:
	double Omega = 1.0,
	       length = 1.1,
		   mass = 7.0,
		   A = 0.15,
		   delta1 = 0.245,
		   delta2 = 0.245;

	double f1( double sinphi, double cosphi ) const;
	double f2( double cosphi ) const;
	double f3( double cospsi ) const;
	double f4( double cosphi ) const;

	double vpDot( double sinphi, double cosphi, double costheta, double cospsi, double vp, double vt ) const;
	double vtDot( double sinphi, double cosphi, double sintheta, double cospsi, double vp, double vt ) const;

	double f1d( double sinphi, double cosphi ) const;
	double f3d( double sinpsi ) const;
};

struct PendulumFamily : Family<Pendulum> {

	PendulumFamily() : Family<Pendulum>(5,1) {}

	Pendulum operator()( const VectorXd& parameter ) const {
		Pendulum P;
		P.Omega = parameter[0];
		return P;
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

typedef RKDynamicalSystem<Pendulum,PendulumWrap> PendulumSolver;

#endif
