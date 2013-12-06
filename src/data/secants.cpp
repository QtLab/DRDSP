#include <DRDSP/data/secants.h>
#include <iostream>
#include <limits>

using namespace std;
using namespace DRDSP;

Secants::Secants() : dimension(0), count(0) {
}

SecantsPreComputed::SecantsPreComputed() : secants(nullptr), weights(nullptr) {
}

SecantsPreComputed::SecantsPreComputed( const SecantsPreComputed& rhs ) : secants(nullptr), weights(nullptr) {
	count = rhs.count;
	dimension = rhs.dimension;
	if( rhs.weights ) {
		weights = new weightType [count];
		for(uint32_t i=0;i<count;i++)
			weights[i] = rhs.weights[i];
	}
	if( rhs.secants ) {
		secants = new VectorXd [count];
		for(uint32_t i=0;i<count;i++)
			secants[i] = rhs.secants[i];
	}
}

SecantsPreComputed::SecantsPreComputed( SecantsPreComputed&& rhs ) {
	count = rhs.count;
	dimension = rhs.dimension;
	weights = rhs.weights;
	secants = rhs.secants;

	rhs.count = 0;
	rhs.dimension = 0;
	rhs.weights = nullptr;
	rhs.secants = nullptr;
}

SecantsPreComputed::~SecantsPreComputed() {
	delete[] secants;
	delete[] weights;
}

SecantsPreComputed& SecantsPreComputed::operator=( const SecantsPreComputed& rhs ) {
	delete[] secants; secants = nullptr;
	delete[] weights; weights = nullptr;
	count = rhs.count;
	dimension = rhs.dimension;
	if( rhs.weights ) {
		weights = new weightType [count];
		for(uint32_t i=0;i<count;i++)
			weights[i] = rhs.weights[i];
	}
	if( rhs.secants ) {
		secants = new VectorXd [count];
		for(uint32_t i=0;i<count;i++)
			secants[i] = rhs.secants[i];
	}
	return *this;
}

SecantsPreComputed& SecantsPreComputed::operator=( SecantsPreComputed&& rhs ) {
	if( this != &rhs ) {
		delete[] secants;
		delete[] weights;
		count = rhs.count;
		dimension = rhs.dimension;
		weights = rhs.weights;
		secants = rhs.secants;

		rhs.count = 0;
		rhs.dimension = 0;
		rhs.weights = nullptr;
		rhs.secants = nullptr;
	}
	return *this;
}

void SecantsData::SetData( const DataSet& dataSet ) {
	data = &dataSet;
	dimension = dataSet.dimension;
	count = ( dataSet.count * (dataSet.count - 1) ) / 2;
}

void SecantsPreComputed::ComputeFromData( const DataSet& data ) {
	dimension = data.dimension;
	count = ( data.count * (data.count - 1) ) / 2;
	size_t bytes = sizeof(double) * data.dimension * count;

	secants = new VectorXd [count];
	uint32_t k = 0;
	for(uint32_t i=0;i<data.count;i++) {
		for(uint32_t j=i+1;j<data.count;j++) {
			secants[k++] = ( data.points[j] - data.points[i] ).normalized();
		}
	}
}

uint32_t SecantsData::GetIndexI( uint32_t k, uint32_t N ) {
	for(uint32_t a=0;a<N;a++) {
		uint32_t Ma = ((a+1)*(2*(N-1)-a))/2 - 1;
		if( Ma >= k ) {
			return a;
		}
	}
	return N-1;
}

uint32_t SecantsData::GetIndexJ( uint32_t k, uint32_t i, uint32_t N ) {
	return (k + 1 + (i*(i-1))/2) - i*(N-2);
}

VectorXd SecantsPreComputed::GetSecant( uint32_t k ) const {
	return secants[k];
}

VectorXd SecantsPreComputed::GetSecantNoNormalize( uint32_t k ) const {
	return secants[k];
}

VectorXd SecantsData::GetSecant( uint32_t k ) const {
	uint32_t i = GetIndexI(k,data->count);
	uint32_t j = GetIndexJ(k,i,data->count);
	return ( data->points[j] - data->points[i] ).normalized();
}

VectorXd SecantsData::GetSecantNoNormalize( uint32_t k ) const {
	uint32_t i = GetIndexI(k,data->count);
	uint32_t j = GetIndexJ(k,i,data->count);
	return data->points[j] - data->points[i];
}

SecantsPreComputed Secants::CullSecants( double tolerance ) const {
	uint32_t culled = 0;
	VectorXd si, sj;
	double dot, test;
	double tolerance2 = tolerance * tolerance;
	
	weightType* cullWeights = new weightType [count];
	for(uint32_t i=0;i<count;i++) {
		cullWeights[i] = 1;
	}

	for(uint32_t i=0;i<count;i++) {
		if( cullWeights[i] == 0 ) continue;
		si = GetSecantNoNormalize(i);
		for(uint32_t j=i+1;j<count;j++) {
			if( cullWeights[j] == 0 ) continue;
			sj = GetSecantNoNormalize(j);
			dot = si.dot(sj);
			test = ( dot * dot ) / ( si.squaredNorm() * sj.squaredNorm() );
			if( test >= tolerance2 ) {
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

	SecantsPreComputed culledSecants;
	culledSecants.count = remain;
	culledSecants.secants = new VectorXd [remain];
	culledSecants.weights = new weightType [remain];
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

SecantsPreComputed SecantsPreComputed::CullSecants( double tolerance ) const {
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

	SecantsPreComputed culledSecants;
	culledSecants.count = remain;
	culledSecants.secants = new VectorXd [remain];
	culledSecants.weights = new weightType [remain];
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

SecantsPreComputed Secants::CullSecantsDegrees( double degrees ) const {
	return CullSecants( cos( degrees * (M_PI/180.0) ) );
}

SecantsPreComputed Secants::CullSecantsRadians( double radians ) const {
	return CullSecants( cos( radians ) );
}
