#include <DRDSP/data/secants.h>
#include <iostream>
#include <limits>

using namespace std;
using namespace DRDSP;

Secants::Secants() : secants(nullptr), count(0), preComputed(false), data(nullptr), weights(nullptr), weighted(false) {
}

Secants::Secants( const Secants& rhs ) : secants(nullptr), count(0), preComputed(false), data(nullptr), weights(nullptr), weighted(false) {
	data = rhs.data;
	count = rhs.count;
	weighted = rhs.weighted;
	if( weighted ) {
		weights = new weightType [count];
		for(uint32_t i=0;i<count;i++)
			weights[i] = rhs.weights[i];
	}
	preComputed = rhs.preComputed;
	if( preComputed ) {
		secants = new VectorXd [count];
		for(uint32_t i=0;i<count;i++)
			secants[i] = rhs.secants[i];
	}
}

Secants::Secants( Secants&& rhs ) {
	delete[] secants;
	delete[] weights;
	data = rhs.data;
	count = rhs.count;
	weighted = rhs.weighted;
	weights = rhs.weights;
	preComputed = rhs.preComputed;
	secants = rhs.secants;

	rhs.data = nullptr;
	rhs.count = 0;
	rhs.weighted = false;
	rhs.weights = nullptr;
	rhs.preComputed = false;
	rhs.secants = nullptr;
}

Secants::~Secants() {
	delete[] secants;
	delete[] weights;
}

void Secants::ComputeFromData( const DataSet& dataSet, size_t preComputeSize ) {
	data = &dataSet;
	dimension = dataSet.dimension;
	uint32_t count = ( dataSet.count * (dataSet.count - 1) ) / 2;
	size_t bytes = sizeof(double) * dataSet.dimension * count;
	if( bytes <= preComputeSize ) {
		PreCompute();
	}
}

uint32_t Secants::GetIndexI( uint32_t k, uint32_t N ) {
	for(uint32_t a=0;a<N;a++) {
		uint32_t Ma = ((a+1)*(2*(N-1)-a))/2 - 1;
		if( Ma >= k ) {
			return a;
		}
	}
}

uint32_t Secants::GetIndexJ( uint32_t k, uint32_t i, uint32_t N ) {
	return (k + 1 + (i*(i-1))/2) - i*(N-2);
}

VectorXd Secants::GetSecant( uint32_t k ) const {
	if( preComputed ) return secants[k];

	uint32_t i = GetIndexI(k,data->count);
	uint32_t j = GetIndexJ(k,i,data->count);

	return ( data->points[j] - data->points[i] ).normalized();
}


void Secants::PreCompute() {
	secants = new VectorXd [count];
	uint32_t k = 0;
	for(uint32_t i=0;i<data->count;i++) {
		for(uint32_t j=i+1;j<data->count;j++) {
			secants[k++] = ( data->points[j] - data->points[i] ).normalized();
		}
	}
}

Secants Secants::CullSecants( double tolerance ) const {
	cout << "Culling Secants..." << endl;
	uint32_t culled = 0;
	VectorXd si;
	double dot;
	double tolerance2 = tolerance * tolerance;
	
	weightType* weights = new weightType [count];
	memset(weights,0,sizeof(weightType)*count);

	for(uint32_t i=0;i<count;i++) {
		if( weights[i] == 0 ) continue;
		si = GetSecant(i);
		for(uint32_t j=i+1;j<count;j++) {
			if( weights[j] == 0 ) continue;
			dot = si.dot(GetSecant(j));
			if( dot*dot >= tolerance2 ) {
				if( weights[i] == numeric_limits<weightType>::max() ) {
					cout << "CullSecants: Overflow" << endl;
				} else {
					weights[i]++;
					weights[j] = 0;
					culled++;
				}
			}
		}
	}
	uint32_t remain = count - culled;

	Secants culledSecants;
	culledSecants.count = remain;
	culledSecants.secants = new VectorXd [remain];
	culledSecants.weights = new weightType [remain];
	culledSecants.preComputed = true;
	culledSecants.dimension = dimension;

	uint32_t j=0;
	for(uint32_t i=0;i<count;i++) {
		if( weights[i] ) {
			culledSecants.secants[j] = GetSecant(i);
			culledSecants.weights[j] = weights[i];
			j++;
		}
	}
	
	cout << "Culled " << culled << " secants. " << remain << " remain ( " << double(100*remain)/count << "% )" << endl;
	return std::move(culledSecants);
}

Secants Secants::CullSecantsDegrees( double degrees ) const {
	return CullSecants( cos( degrees * (M_PI/180.0) ) );
}

