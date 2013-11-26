#ifndef INCLUDED_GOODFELLOW
#define INCLUDED_GOODFELLOW
#include <DRDSP/dynamics/model.h>

using namespace DRDSP;

struct Goodfellow : ModelParameterized {
	MatrixXd adjacency;
	VectorXd mu;
	double omega, a, b, c, d;
	uint32_t N;

	Goodfellow();
	explicit Goodfellow( uint32_t numCompartments );
	void Create( uint32_t numCompartments );
	VectorXd VectorField( const VectorXd& state, const VectorXd& parameter );
	MatrixXd Partials( const VectorXd& state, const VectorXd& parameter );

	
	double poly( const VectorXd& state, uint32_t i ) {
		return mu(i) - a * r2(state,i) + b * r4(state,i) - c * r6(state,i);
	}

	double dpolydx( const VectorXd& state, uint32_t i, uint32_t j ) {
		return - a * dr2dx(state,i,j) + b * dr4dx(state,i,j) - c * dr6dx(state,i,j);
	}

	double dpolydy( const VectorXd& state, uint32_t i, uint32_t j ) {
		return - a * dr2dy(state,i,j) + b * dr4dy(state,i,j) - c * dr6dy(state,i,j);
	}

	double x( const VectorXd& state, uint32_t i ) {
		return state(i);
	}

	double y( const VectorXd& state, uint32_t i ) {
		return state(N+i);
	}

	double r2( const VectorXd& state, uint32_t i ) {
		return x(state,i)*x(state,i) + y(state,i)*y(state,i);
	}

	double r4( const VectorXd& state, uint32_t i ) {
		return r2(state,i)*r2(state,i);
	}

	double r6( const VectorXd& state, uint32_t i ) {
		return r4(state,i)*r2(state,i);
	}

	double dr2dx( const VectorXd& state, uint32_t i, uint32_t j ) {
		return delta(i,j) * 2.0 * x(state,i);
	}

	double dr4dx( const VectorXd& state, uint32_t i, uint32_t j ) {
		return delta(i,j) * 4.0 * x(state,i) * r2(state,i);
	}

	double dr6dx( const VectorXd& state, uint32_t i, uint32_t j ) {
		return delta(i,j) * 6.0 * x(state,i) * r4(state,i);
	}

	double dr2dy( const VectorXd& state, uint32_t i, uint32_t j ) {
		return delta(i,j) * 2.0 * y(state,i);
	}

	double dr4dy( const VectorXd& state, uint32_t i, uint32_t j ) {
		return delta(i,j) * 4.0 * y(state,i) * r2(state,i);
	}

	double dr6dy( const VectorXd& state, uint32_t i, uint32_t j ) {
		return delta(i,j) * 6.0 * y(state,i) * r4(state,i);
	}

	double delta( uint32_t i, uint32_t j ) {
		return (i==j)?1.0:0.0;
	}

};

#endif
