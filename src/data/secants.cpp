#include <DRDSP/data/secants.h>
#include <iostream>
#include <limits>

using namespace std;
using namespace DRDSP;

Secants::Secants() : dimension(0), count(0) {}

void SecantsData::SetData( const DataSet& dataSet ) {
	data = &dataSet;
	dimension = dataSet.dimension;
	size_t N = dataSet.points.size();
	count = ( N * (N - 1) ) / 2;
}

void SecantsPreComputed::ComputeFromData( const DataSet& data ) {
	dimension = data.dimension;
	size_t N = data.points.size();
	count = ( N * (N - 1) ) / 2;
	size_t bytes = sizeof(double) * data.dimension * count;

	secants.resize(count);
	size_t k = 0;
	for(size_t i=0;i<N;i++) {
		for(size_t j=i+1;j<N;j++) {
			secants[k++] = ( data.points[j] - data.points[i] ).normalized();
		}
	}
}

size_t SecantsData::GetIndexI( size_t k, size_t N ) {
	for(size_t a=0;a<N;a++) {
		size_t Ma = ((a+1)*(2*(N-1)-a))/2 - 1;
		if( Ma >= k ) {
			return a;
		}
	}
	return N-1;
}

size_t SecantsData::GetIndexJ( size_t k, size_t i, size_t N ) {
	return (k + 1 + (i*(i-1))/2) - i*(N-2);
}

VectorXd SecantsPreComputed::GetSecant( size_t k ) const {
	return secants[k];
}

VectorXd SecantsPreComputed::GetSecantNoNormalize( size_t k ) const {
	return secants[k];
}

VectorXd SecantsData::GetSecant( size_t k ) const {
	size_t i = GetIndexI(k,data->points.size());
	size_t j = GetIndexJ(k,i,data->points.size());
	return ( data->points[j] - data->points[i] ).normalized();
}

VectorXd SecantsData::GetSecantNoNormalize( size_t k ) const {
	size_t i = GetIndexI(k,data->points.size());
	size_t j = GetIndexJ(k,i,data->points.size());
	return data->points[j] - data->points[i];
}

SecantsPreComputed Secants::CullSecants( double tolerance ) const {
	size_t culled = 0;
	VectorXd si, sj;
	double dot, test;
	double tolerance2 = tolerance * tolerance;
	
	vector<weightType> cullWeights(count);
	for(size_t i=0;i<count;i++) {
		cullWeights[i] = 1;
	}

	for(size_t i=0;i<count;i++) {
		if( cullWeights[i] == 0 ) continue;
		si = GetSecantNoNormalize(i);
		for(size_t j=i+1;j<count;j++) {
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
	size_t remain = count - culled;

	SecantsPreComputed culledSecants;
	culledSecants.count = remain;
	culledSecants.secants.resize(remain);
	culledSecants.weights.resize(remain);
	culledSecants.dimension = dimension;

	size_t j=0;
	for(size_t i=0;i<count;i++) {
		if( cullWeights[i] ) {
			culledSecants.secants[j] = GetSecant(i);
			culledSecants.weights[j] = cullWeights[i];
			j++;
		}
	}

	cout << "Culled " << culled << " secants. " << remain << " remain ( " << double(100*remain)/count << "% )" << endl;
	return culledSecants;
}

SecantsPreComputed SecantsPreComputed::CullSecants( double tolerance ) const {
	size_t culled = 0;
	VectorXd si;
	double dot;
	double tolerance2 = tolerance * tolerance;
	
	vector<weightType> cullWeights(count);
	for(uint32_t i=0;i<count;i++) {
		cullWeights[i] = 1;
	}

	for(size_t i=0;i<count;i++) {
		if( cullWeights[i] == 0 ) continue;
		si = GetSecant(i);
		for(size_t j=i+1;j<count;j++) {
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
	size_t remain = count - culled;

	SecantsPreComputed culledSecants;
	culledSecants.count = remain;
	culledSecants.secants.resize(remain);
	culledSecants.weights.resize(remain);
	culledSecants.dimension = dimension;

	size_t j=0;
	for(size_t i=0;i<count;i++) {
		if( cullWeights[i] ) {
			culledSecants.secants[j] = GetSecant(i);
			culledSecants.weights[j] = cullWeights[i];
			j++;
		}
	}

	cout << "Culled " << culled << " secants. " << remain << " remain ( " << double(100*remain)/count << "% )" << endl;
	return culledSecants;
}

SecantsPreComputed Secants::CullSecantsDegrees( double degrees ) const {
	return CullSecants( cos( degrees * (M_PI/180.0) ) );
}

SecantsPreComputed Secants::CullSecantsRadians( double radians ) const {
	return CullSecants( cos( radians ) );
}
