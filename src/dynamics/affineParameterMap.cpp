#include <DRDSP/dynamics/affineParameterMap.h>

using namespace DRDSP;

AffineParameterMap::AffineParameterMap( uint32_t dim, uint32_t numRBFs, uint32_t paramDim ) :
	N( dim + numRBFs ),
	R( N * (uint32_t)(paramDim+1) ),
	dimension(dim),
	parameterDimension(paramDim)
{
	coeffs.setZero(dimension,R);
}

MatrixXd AffineParameterMap::operator()( const VectorXd &parameter ) const {
	return coeffs * GetLambda(parameter);
}

MatrixXd AffineParameterMap::GetLambda( const VectorXd& parameter ) const {
	MatrixXd Lambda;
	Lambda.setZero(R,N);
	Lambda.block(0,0,N,N).setIdentity();
	for(uint32_t i=0;i<parameterDimension;i++) {
		Lambda.block(N*(i+1),0,N,N).setIdentity();
		Lambda.block(N*(i+1),0,N,N) *= parameter(i);
	}
	return Lambda;
}
