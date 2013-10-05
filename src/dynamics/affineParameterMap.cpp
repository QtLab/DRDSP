#include <DRDSP/dynamics/affineParameterMap.h>

using namespace DRDSP;

void AffineParameterMap::Init( uint32_t dim, uint32_t numRBFs, uint32_t paramDim ) {
	dimension = dim;
	N = dimension + numRBFs;
	parameterDimension = paramDim;
	R = N * (parameterDimension+1);
	coeffs.setZero(dimension,R);
}

MatrixXd AffineParameterMap::Evaluate( const VectorXd &x ) const {
	return coeffs * GetLambda(x);
}

MatrixXd AffineParameterMap::GetLambda( const VectorXd &x ) const {
	MatrixXd Lambda;
	Lambda.setZero(R,N);
	Lambda.block(0,0,N,N).setIdentity();
	for(uint32_t i=0;i<parameterDimension;i++) {
		Lambda.block(N*(i+1),0,N,N).setIdentity();
		Lambda.block(N*(i+1),0,N,N) *= x(i);
	}
	return Lambda;
}
