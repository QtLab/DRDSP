#include <DRDSP/data/secants.h>
#include <DRDSP/geometry/trig.h>
#include <iostream>
#include <limits>

using namespace std;
using namespace DRDSP;

VectorXd Secants::GetSecant( size_t k ) const {
	return secants[k];
}

Secants Secants::CullSecantsAngle( Degreesd degrees ) const {
	return CullSecants( cos( Radiansd(degrees) ) );
}

Secants Secants::CullSecantsAngle( Radiansd radians ) const {
	return CullSecants( cos( radians ) );
}

Secants& Secants::ComputeFromData( const DataSet& data ) {
	dimension = data.dimension;
	size_t N = data.points.size();
	count = ( N * (N - 1) ) / 2;

	secants.resize(count);
	VectorXd diff;
	size_t k = 0;
	for(size_t i=0;i<N;++i) {
		for(size_t j=i+1;j<N;++j) {
			diff = data.points[j] - data.points[i];
			double norm2 = diff.squaredNorm();
			if( norm2 != 0.0 ) {
				secants[k++] = diff / sqrt(norm2);
			}
		}
	}
	count = k;
	secants.resize(count);
	return *this;
}

Secants& Secants::ComputeFromData( const DataSet& data, double tolerance ) {
	size_t culled = 0;
	VectorXd si, sj;
	double dot;
	double tolerance2 = tolerance * tolerance;
	size_t N = data.points.size();
	size_t M = ( N * ( N - 1 ) ) / 2;

	vector<weightType> cullWeights(M);
	for(uint32_t i=0;i<M;++i) {
		cullWeights[i] = 1;
	}
	size_t ia = 0, ib = 1, ja, jb;
	for(size_t i=0;i<M;++i) {
		if( cullWeights[i] == 0 ) {
			++ib;
			if( ib == N ) {
				++ia;
				ib = ia + 1;
			}
			continue;
		}
		si = data.points[ib] - data.points[ia];
		ja = ia;
		jb = ib;
		++jb;
		if( jb == N ) {
			++ja;
			jb = ja + 1;
		}
		for(size_t j=i+1;j<M;++j) {
			if( cullWeights[j] == 0 ) {
				++jb;
				if( jb == N ) {
					++ja;
					jb = ja + 1;
				}
				continue;
			}
			sj = data.points[jb] - data.points[ja];
			dot = si.dot(sj);
			if( dot*dot >= tolerance2 * si.squaredNorm() * sj.squaredNorm() ) {
				if( cullWeights[i] == numeric_limits<weightType>::max() ) {
					cout << "CullSecants: Overflow" << endl;
				} else {
					cullWeights[i]++;
					cullWeights[j] = 0;
					culled++;
				}
			}
			++jb;
			if( jb == N ) {
				++ja;
				jb = ja + 1;
			}
		}
		++ib;
		if( ib == N ) {
			++ia;
			ib = ia + 1;
		}
	}
	size_t remain = M - culled;

	count = remain;
	secants.resize(remain);
	weights.resize(remain);
	dimension = data.dimension;

	size_t j=0;
	ia = 0; ib = 1;
	for(size_t i=0;i<M;++i) {
		if( cullWeights[i] ) {
			secants[j] = (data.points[ib] - data.points[ia]).normalized();
			weights[j] = cullWeights[i];
			++j;
		}
		++ib;
		if( ib == N ) {
			++ia;
			ib = ia + 1;
		}
	}

	cout << "Culled " << culled << " secants. " << remain << " remain ( " << double(100*remain)/M << "% )" << endl;
	return *this;
}

Secants Secants::CullSecants( double tolerance ) const {
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

	Secants culledSecants;
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


vector<Secants> DRDSP::ComputeSecants( const DataSystem& data ) {
	vector<Secants> secants( data.numParameters );
	for(uint32_t i=0;i<data.numParameters;++i)
		secants[i].ComputeFromData( data.dataSets[i] );
	return secants;
}

vector<Secants> DRDSP::ComputeSecants( const DataSystem& data, uint32_t numThreads ) {
	vector<Secants> secants( data.numParameters );
	vector<future<void>> futures(numThreads);
			
	for(uint32_t i=0;i<data.numParameters;i+=numThreads) {
		uint32_t N = min( data.numParameters - i, numThreads );
		for(uint32_t j=0;j<N;++j) {
			futures[j] = async( launch::async,
				[]( Secants& secants, const DataSet& dataSet ) {
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

vector<Secants> DRDSP::CullSecants( const vector<Secants>& secants, Degreesd degrees ) {
	vector<Secants> culledSecants( secants.size() );
	double tolerance = cos( Radiansd(degrees) );
	for(size_t i=0;i<secants.size();++i)
		culledSecants[i] = secants[i].CullSecants( tolerance );
	return culledSecants;
}

vector<Secants> DRDSP::CullSecants( const vector<Secants>& secants, Degreesd degrees, uint32_t numThreads ) {
	vector<Secants> culledSecants( secants.size() );
	vector<future<void>> futures(numThreads);
	uint32_t numParams = (uint32_t)secants.size();
	double tolerance = cos( Radiansd(degrees) );
	for(uint32_t i=0;i<numParams;i+=numThreads) {
		uint32_t N = min( numParams - i, numThreads );
		for(uint32_t j=0;j<N;++j) {
			futures[j] = async( launch::async,
				[tolerance]( Secants& culled, const Secants& sec ) {
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

vector<Secants> DRDSP::ComputeSecants( const DataSystem& data, Degreesd degrees, uint32_t numThreads ) {
	vector<Secants> culledSecants( data.numParameters );
	vector<future<void>> futures(numThreads);
	uint32_t numParams = data.numParameters;
	double tolerance = cos( Radiansd(degrees) );
	for(uint32_t i=0;i<numParams;i+=numThreads) {
		uint32_t N = min( numParams - i, numThreads );
		for(uint32_t j=0;j<N;++j) {
			futures[j] = async( launch::async,
				[tolerance]( Secants& culled, const DataSet& dataSet ) {
					culled.ComputeFromData( dataSet, tolerance );
				},
				ref(culledSecants[i+j]), cref(data.dataSets[i+j])
			);
		}
		for(uint32_t j=0;j<N;++j) {
			futures[j].wait();
		}
	}
	return culledSecants;
}
