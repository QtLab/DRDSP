#ifndef INCLUDED_DYNAMICS_MODEL_RBF
#define INCLUDED_DYNAMICS_MODEL_RBF
#include "../types.h"
#include "model.h"
#include "radial_basis.h"
#include "../data/aabb.h"

namespace DRDSP {

	struct ModelRBF : Model {
		MatrixXd linear;
		VectorXd* weights;
		RadialFunction* rbfs;
		uint16_t numRBFs;

		ModelRBF();
		ModelRBF( const ModelRBF& rhs );
		ModelRBF( ModelRBF&& rhs );
		ModelRBF( uint16_t dim, uint16_t nRBFs );
		ModelRBF& operator=( const ModelRBF& rhs );
		ModelRBF& operator=( ModelRBF&& rhs );
		void Create( uint16_t dim, uint16_t nRBFs );
		void Destroy();
		VectorXd VectorField( const VectorXd& x );
		MatrixXd Partials( const VectorXd &x );
		void SetCentresRandom( const AABB& box );
		void SetRBFType( const Function& f );
		void LoadCentresText( const char* filename );
		void LoadCentresBinary( const char* filename );
		void OutputText( const char *filename ) const;
	};

}

#endif
