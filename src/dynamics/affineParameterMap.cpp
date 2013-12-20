#include <DRDSP/dynamics/affineParameterMap.h>

using namespace DRDSP;

void AffineParameterMap::Init( uint16_t dim, uint16_t numRBFs, uint8_t paramDim ) {
	dimension = dim;
	N = dimension + numRBFs;
	parameterDimension = paramDim;
	R = N * (uint32_t)(parameterDimension+1);
	coeffs.setZero(dimension,R);
}

MatrixXd AffineParameterMap::Evaluate( const VectorXd &parameter ) const {
	return coeffs * GetLambda(parameter);
}

MatrixXd AffineParameterMap::GetLambda( const VectorXd &parameter ) const {
	MatrixXd Lambda;
	Lambda.setZero(R,N);
	Lambda.block(0,0,N,N).setIdentity();
	for(uint8_t i=0;i<parameterDimension;i++) {
		Lambda.block(N*(i+1),0,N,N).setIdentity();
		Lambda.block(N*(i+1),0,N,N) *= parameter(i);
	}
	return Lambda;
}
