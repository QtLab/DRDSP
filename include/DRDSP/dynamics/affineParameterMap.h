#ifndef INCLUDED_DYNAMICS_AFFINEPARAMETERMAP
#define INCLUDED_DYNAMICS_AFFINEPARAMETERMAP
#include "../types.h"

namespace DRDSP {

	/*!
	 * \brief An affine map from the original parameter space to the reduced parameter space
	 */
	struct AffineParameterMap {
		uint32_t R,                 //!< Width of the coefficient matrix
			     N;                 //!< Width of the matrix version of the parameter returned by GetLambda()
		uint16_t dimension;         //!< Dimension of the reduced state space
		uint8_t parameterDimension; //!< Dimension of the original parameter space
		MatrixXd coeffs;            //!< Coefficients describing the affine map
	
		void Init( uint16_t dim, uint16_t numRBFs, uint8_t paramDim );
		MatrixXd Evaluate( const VectorXd &parameter ) const;          //!< Evaluate the parameter map
		MatrixXd GetLambda( const VectorXd &parameter ) const;         //!< Get the matrix version of the original parameter
	};

}

#endif

