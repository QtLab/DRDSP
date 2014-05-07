#ifndef INCLUDED_DYNAMICS_AFFINEPARAMETERMAP
#define INCLUDED_DYNAMICS_AFFINEPARAMETERMAP
#include "../types.h"

namespace DRDSP {

	/*!
	 * \brief An affine map from the original parameter space to the reduced parameter space
	 */
	struct AffineParameterMap {
		uint32_t N,                 //!< Width of the matrix version of the parameter returned by GetLambda()
			     R;                 //!< Width of the coefficient matrix
		uint32_t dimension;         //!< Dimension of the reduced state space
		uint32_t parameterDimension; //!< Dimension of the original parameter space
		MatrixXd coeffs;            //!< Coefficients describing the affine map
	
		AffineParameterMap() = default;
		AffineParameterMap( uint32_t dim, uint32_t numRBFs, uint32_t paramDim );
		MatrixXd operator()( const VectorXd &parameter ) const;        //!< Evaluate the parameter map
		MatrixXd GetLambda( const VectorXd &parameter ) const;         //!< Get the matrix version of the original parameter
	};

}

#endif

