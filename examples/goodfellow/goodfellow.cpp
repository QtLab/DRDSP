#include "goodfellow.h"

Goodfellow::Goodfellow() : Goodfellow(100) {}

Goodfellow::Goodfellow( uint32_t numCompartments ) :
	Model<>(2*numCompartments),
	N(numCompartments),
	omega(20.0),
	a(2.0),
	b(3.0/2.0),
	c(1.0/3.0),
	d(10.0)
{
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
		for(uint32_t i=2;i<N;++i)
			mu(i) = muMax - ((muMax-muMin)*i)/(N-1);
	}

	adjacency.setZero(N,N);
	for(uint32_t i=0;i<N;++i)
		for(uint32_t j=0;j<N;++j)
			adjacency(i,j) = (i==j)?0.0:1.0;
}

