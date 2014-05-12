#include <DRDSP/data/secants.h>
#include <DRDSP/geometry/trig.h>
#include <iostream>
#include <limits>

using namespace std;
using namespace DRDSP;

Secants::Secants() : dimension(0), count(0) {}

SecantsPreComputed Secants::CullSecantsAngle( Degreesd degrees ) const {
	return CullSecants( cos( Radiansd(degrees) ) );
}

SecantsPreComputed Secants::CullSecantsAngle( Radiansd radians ) const {
	return CullSecants( cos( radians ) );
}

SecantsPreComputed Secants::CullSecants( double tolerance ) const {
	size_t culled = 0;
	VectorXd si, sj;
	double dot, test;
	double tolerance2 = tolerance * tolerance;
	
	vector<weightType> cullWeights(count);
	for(size_t i=0;i<count;++i) {
		cullWeights[i] = 1;
	}

	for(size_t i=0;i<count;++i) {
		if( cullWeights[i] == 0 ) continue;
		si = GetSecantNoNormalize(i);
		for(size_t j=i+1;j<count;++j) {
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
	for(size_t i=0;i<count;++i) {
		if( cullWeights[i] ) {
			culledSecants.secants[j] = GetSecant(i);
			culledSecants.weights[j] = cullWeights[i];
			++j;
		}
	}

	cout << "Culled " << culled << " secants. " << remain << " remain ( " << double(100*remain)/count << "% )" << endl;
	return culledSecants;
}


void SecantsPreComputed::ComputeFromData( const DataSet& data ) {
	dimension = data.dimension;
	size_t N = data.points.size();
	count = ( N * (N - 1) ) / 2;

	secants.resize(count);
	size_t k = 0;
	for(size_t i=0;i<N;++i) {
		for(size_t j=i+1;j<N;++j) {
			secants[k++] = ( data.points[j] - data.points[i] ).normalized();
		}
	}
}

VectorXd SecantsPreComputed::GetSecant( size_t k ) const {
	return secants[k];
}

VectorXd SecantsPreComputed::GetSecantNoNormalize( size_t k ) const {
	return secants[k];
}

SecantsPreComputed SecantsPreComputed::CullSecants( double tolerance ) const {
	size_t culled = 0;
	VectorXd si;
	double dot;
	double tolerance2 = tolerance * tolerance;
	
	vector<weightType> cullWeights(count);
	for(uint32_t i=0;i<count;++i) {
		cullWeights[i] = 1;
	}

	for(size_t i=0;i<count;++i) {
		if( cullWeights[i] == 0 ) continue;
		si = GetSecant(i);
		for(size_t j=i+1;j<count;++j) {
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
	for(size_t i=0;i<count;++i) {
		if( cullWeights[i] ) {
			culledSecants.secants[j] = GetSecant(i);
			culledSecants.weights[j] = cullWeights[i];
			++j;
		}
	}

	cout << "Culled " << culled << " secants. " << remain << " remain ( " << double(100*remain)/count << "% )" << endl;
	return culledSecants;
}

void SecantsData::SetData( const DataSet& dataSet ) {
	data = &dataSet;
	dimension = dataSet.dimension;
	size_t N = dataSet.points.size();
	count = ( N * (N - 1) ) / 2;
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

size_t SecantsData::GetIndexI( size_t k, size_t N ) {
	for(size_t a=0;a<N;++a) {
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

vector<SecantsPreComputed> DRDSP::ComputeSecants( const DataSystem& data ) {
	vector<SecantsPreComputed> secants( data.numParameters );
	for(uint32_t i=0;i<data.numParameters;++i)
		secants[i].ComputeFromData( data.dataSets[i] );
	return secants;
}

vector<SecantsPreComputed> DRDSP::ComputeSecants( const DataSystem& data, uint32_t numThreads ) {
	vector<SecantsPreComputed> secants( data.numParameters );
	vector<future<void>> futures(numThreads);
			
	for(uint32_t i=0;i<data.numParameters;i+=numThreads) {
		uint32_t N = min( data.numParameters - i, numThreads );
		for(uint32_t j=0;j<N;++j) {
			futures[j] = async( launch::async,
				[]( SecantsPreComputed& secants, const DataSet& dataSet ) {
					secants.ComputeFromData( dataSet );
				},
				ref(secants[i+j]), cref(data.dataSets[i+j])
			);
		}
		for(uint32_t j=0;j<N;++j) {
			futures[j].wait();
		}
	}
	return secants;
}

vector<SecantsPreComputed> DRDSP::CullSecants( const vector<SecantsPreComputed>& secants, Degreesd degrees ) {
	vector<SecantsPreComputed> culledSecants( secants.size() );
	double tolerance = cos( Radiansd(degrees) );
	for(size_t i=0;i<secants.size();++i)
		culledSecants[i] = secants[i].CullSecants( tolerance );
	return culledSecants;
}

vector<SecantsPreComputed> DRDSP::CullSecants( const vector<SecantsPreComputed>& secants, Degreesd degrees, uint32_t numThreads ) {
	vector<SecantsPreComputed> culledSecants( secants.size() );
	vector<future<void>> futures(numThreads);
	uint32_t numParams = (uint32_t)secants.size();
	double tolerance = cos( Radiansd(degrees) );
	for(uint32_t i=0;i<numParams;i+=numThreads) {
		uint32_t N = min( numParams - i, numThreads );
		for(uint32_t j=0;j<N;++j) {
			futures[j] = async( launch::async,
				[tolerance]( SecantsPreComputed& culled, const SecantsPreComputed& sec ) {
					culled = sec.CullSecants( tolerance );
				},
				ref(culledSecants[i+j]), cref(secants[i+j])
			);
		}
		for(uint32_t j=0;j<N;++j) {
			futures[j].wait();
		}
	}
	return culledSecants;
}
