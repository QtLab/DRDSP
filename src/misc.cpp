#include <DRDSP/misc.h>

using namespace std;

vector<VectorXd> DRDSP::ParameterList( const VectorXd& pMin, const VectorXd& pMax, uint32_t N ) {
	vector<VectorXd> v(N);
	VectorXd pDelta = ( pMax - pMin ) / ( N - 1 );
	v[0] = pMin;
	for(uint32_t i=1;i<N;++i) {
		v[i] = v[i-1] + pDelta;
	}
	return v;
}

vector<VectorXd> DRDSP::ParameterList( double pMin, double pMax, uint32_t N ) {
	vector<VectorXd> v(N);
	double pDelta = ( pMax - pMin ) / ( N - 1 );
	VectorXd temp(1);
	temp[0] = pMin;
	v[0] = temp;
	for(uint32_t i=1;i<N;++i) {
		temp[0] = pDelta;
		v[i] = v[i-1] + temp;
	}
	return v;
}

void DRDSP::SetPointsRandom( vector<VectorXd>& points, const AABB& box, mt19937& mt ) {
	static const uniform_real_distribution<double> dist;
	VectorXd diff = box.bMax - box.bMin;
	for( auto& x : points ) {
		for(int64_t j=0;j<x.size();++j) {
			x[j] = box.bMin(j) + diff(j) * dist(mt);
		}
	}
}

