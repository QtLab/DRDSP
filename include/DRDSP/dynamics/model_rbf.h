#ifndef INCLUDED_DYNAMICS_MODEL_RBF
#define INCLUDED_DYNAMICS_MODEL_RBF
#include "../types.h"
#include "radial_basis.h"
#include "../data/aabb.h"

namespace DRDSP {

	struct ModelRBF {
		MatrixXd linear;
		VectorXd* weights;
		RadialFunction* rbfs;
		uint16_t dimension, numRBFs;

		ModelRBF();
		ModelRBF( uint16_t dim, uint16_t nRBFs );
		void Create( uint16_t dim, uint16_t nRBFs );
		void Destroy();
		VectorXd VectorField( const VectorXd& x ) const;
		MatrixXd VectorFieldDerivative( const VectorXd &x ) const;
		void SetCentresRandom( const AABB& box );
		void SetRBFType( const Function& f );
		void LoadCentresText( const char* filename );
		void LoadCentresBinary( const char* filename );
		void OutputText( const char *filename ) const;
	};

}

#endif
