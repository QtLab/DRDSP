#ifndef INCLUDED_DYNAMICS_AFFINEPARAMETERMAP
#define INCLUDED_DYNAMICS_AFFINEPARAMETERMAP
#include "../types.h"

namespace DRDSP {

	struct AffineParameterMap {
		uint32_t R, N;
		uint16_t dimension;
		uint8_t parameterDimension;
		MatrixXd coeffs;
	
		void Init( uint16_t dim, uint16_t numRBFs, uint8_t paramDim );
		MatrixXd Evaluate( const VectorXd &x ) const;
		MatrixXd GetLambda( const VectorXd &x ) const;
	};

}

#endif

