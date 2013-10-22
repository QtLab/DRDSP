#include <DRDSP/data/secants.h>
#include <iostream>
#include <limits>

using namespace std;
using namespace DRDSP;

Secants::Secants() : secants(nullptr), dimension(0), count(0), preComputed(false), data(nullptr), weights(nullptr), weighted(false) {
}

Secants::Secants( const Secants& rhs ) : secants(nullptr), weights(nullptr) {
	data = rhs.data;
	count = rhs.count;
	dimension = rhs.dimension;
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
	data = rhs.data;
	count = rhs.count;
	dimension = rhs.dimension;
	weighted = rhs.weighted;
	weights = rhs.weights;
	preComputed = rhs.preComputed;
	secants = rhs.secants;

	rhs.data = nullptr;
	rhs.count = 0;
	rhs.dimension = 0;
	rhs.weighted = false;
	rhs.weights = nullptr;
	rhs.preComputed = false;
	rhs.secants = nullptr;
}

Secants::~Secants() {
	delete[] secants;
	delete[] weights;
}

Secants& Secants::operator=( const Secants& rhs ) {
	delete[] secants; secants = nullptr;
	delete[] weights; weights = nullptr;
	count = rhs.count;
	data = rhs.data;
	dimension = rhs.dimension;
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
	return *this;
}

Secants& Secants::operator=( Secants&& rhs ) {
	if( this != &rhs ) {
		delete[] secants;
		delete[] weights;
		data = rhs.data;
		count = rhs.count;
		dimension = rhs.dimension;
		weighted = rhs.weighted;
		weights = rhs.weights;
		preComputed = rhs.preComputed;
		secants = rhs.secants;

		rhs.data = nullptr;
		rhs.count = 0;
		rhs.dimension = 0;
		rhs.weighted = false;
		rhs.weights = nullptr;
		rhs.preComputed = false;
		rhs.secants = nullptr;
	}
	return *this;
}

void Secants::ComputeFromData( const DataSet& dataSet, size_t preComputeSize ) {
	data = &dataSet;
	dimension = dataSet.dimension;
	count = ( dataSet.count * (dataSet.count - 1) ) / 2;
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
	return N-1;
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
	cout << "Precomputing Secants..." << endl;
	secants = new VectorXd [count];
	uint32_t k = 0;
	for(uint32_t i=0;i<data->count;i++) {
		for(uint32_t j=i+1;j<data->count;j++) {
			secants[k++] = ( data->points[j] - data->points[i] ).normalized();
		}
	}
	preComputed = true;
}

Secants Secants::CullSecants( double tolerance ) const {
	cout << "Culling Secants..." << endl;
	uint32_t culled = 0;
	VectorXd si;
	double dot;
	double tolerance2 = tolerance * tolerance;
	
	weightType* cullWeights = new weightType [count];
	for(uint32_t i=0;i<count;i++) {
		cullWeights[i] = 1;
	}

	for(uint32_t i=0;i<count;i++) {
		if( cullWeights[i] == 0 ) continue;
		si = GetSecant(i);
		for(uint32_t j=i+1;j<count;j++) {
			if( cullWeights[j] == 0 ) continue;
			dot = si.dot(GetSecant(j));
			if( dot*dot >= tolerance2 ) {
				if( cullWeights[i] == numeric_limits<weightType>::max() ) {
					cout << "CullSecants: Overflow" << endl;
				} else {
					cullWeights[i]++;
					cullWeights[j] = 0;
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
		if( cullWeights[i] ) {
			culledSecants.secants[j] = GetSecant(i);
			culledSecants.weights[j] = cullWeights[i];
			j++;
		}
	}
	
	delete[] cullWeights;

	cout << "Culled " << culled << " secants. " << remain << " remain ( " << double(100*remain)/count << "% )" << endl;
	return std::move(culledSecants);
}

Secants Secants::CullSecantsDegrees( double degrees ) const {
	return CullSecants( cos( degrees * (M_PI/180.0) ) );
}

