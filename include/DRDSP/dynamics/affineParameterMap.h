#ifndef INCLUDED_DYNAMICS_AFFINEPARAMETERMAP
#define INCLUDED_DYNAMICS_AFFINEPARAMETERMAP
#include "../types.h"

namespace DRDSP {

	struct AffineParameterMap {
		uint32_t dimension, R, N, parameterDimension;
		MatrixXd coeffs;
	
		void Init( uint32_t dim, uint32_t numRBFs, uint32_t paramDim );
		MatrixXd Evaluate( const VectorXd &x ) const;
		MatrixXd GetLambda( const VectorXd &x ) const;
	};

}

#endif

