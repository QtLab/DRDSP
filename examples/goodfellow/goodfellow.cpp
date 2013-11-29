#include "goodfellow.h"

Goodfellow::Goodfellow() : N(0), omega(1.0), a(2.0), b(3.0/2.0), c(1.0/3.0), d(1.0) {}

Goodfellow::Goodfellow( uint32_t numCompartments ) : ModelParameterized(2*numCompartments,1), N(numCompartments), omega(1.0), a(2.0), b(3.0/2.0), c(1.0/3.0), d(1.0) {
	Create(N);
}

void Goodfellow::Create( uint32_t numCompartments ) {
	N = numCompartments;
	dimension = 2*N;

	mu.setZero(N);
	double muMin = 0.2;
	double muMax = 0.3;
	mu(0) = 0.6;
	if( N >= 2 ) {
		mu(1) = muMax;
		for(uint32_t i=2;i<N;i++)
			mu(i) = muMax - ((muMax-muMin)*i)/(N-1);
	}

	adjacency.setZero(N,N);
	for(uint32_t i=0;i<N;i++)
		for(uint32_t j=0;j<N;j++)
			adjacency(i,j) = (i==j)?0.0:1.0;
}

VectorXd Goodfellow::VectorField( const VectorXd& state, const VectorXd& parameter ) {
	VectorXd res(2*N);
	for(uint32_t i=0;i<N;i++) {
		double temp1 = omega - d * r2(state,i);
		double temp2 = poly(state,i);
		res(i) = y(state,i) * temp1 + x(state,i) * temp2;
		res(N+i) = -x(state,i) * temp1 + y(state,i) * temp2;
	}
	res.head(N) += parameter(0) * (adjacency * state.head(N)) / N;
	return std::move(res);
}

MatrixXd Goodfellow::Partials( const VectorXd& state, const VectorXd& parameter ) {
	MatrixXd res(2*N,2*N);
	for(uint32_t i=0;i<N;i++) {
		for(uint32_t j=0;j<N;j++) {
			res(i,j) = y(state,i) * (-d * dr2dx(state,i,j)) + delta(i,j) * poly(state,i) + x(state,i) * dpolydx(state,i,j) + parameter(0) * adjacency(i,j) / N;
			res(N+i,j) = -delta(i,j) * ( omega - d * r2(state,i) ) -x(state,i) * (-d * dr2dx(state,i,j)) + y(state,i) * dpolydx(state,i,j);
			res(i,N+j) = delta(i,j) * ( omega - d * r2(state,i) ) + y(state,i) * (-d * dr2dy(state,i,j)) + x(state,i) * dpolydy(state,i,j);
			res(N+i,N+j) = -x(state,i) * (-d * dr2dy(state,i,j)) + delta(i,j) * poly(state,i) + y(state,i) * dpolydy(state,i,j);
		}
	}
	return std::move(res);
}

