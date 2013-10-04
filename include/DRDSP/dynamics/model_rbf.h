#ifndef INCLUDED_DYNAMICS_MODEL_RBF
#define INCLUDED_DYNAMICS_MODEL_RBF
#include "../types.h"
#include "radial_basis.h"

namespace DRDSP {

	struct ModelRBF {
		MatrixXd linear;
		VectorXd* weights;
		RadialFunction* rbfs;
		uint32_t dimension;
		uint16_t numRBFs;

		ModelRBF();
		ModelRBF( uint32_t dim, uint32_t nRBFs );
		void Create( uint32_t dim, uint32_t nRBFs );
		void Destroy();
		VectorXd VectorField( const VectorXd& x ) const;
		MatrixXd VectorFieldDerivative( const VectorXd &x ) const;
		void SetCentresRandom( const VectorXd& minBounds, const VectorXd& maxBounds );
		void SetRBFType( const Function& f );
	};

}

#endif
