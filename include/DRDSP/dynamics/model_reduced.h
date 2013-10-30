#ifndef INCLUDED_MODEL_REDUCED
#define INCLUDED_MODEL_REDUCED
#include <fstream>
#include "model.h"
#include "../types.h"
#include "affineParameterMap.h"
#include "model_rbf.h"
#include "reduced_data_system.h"

using namespace std;

namespace DRDSP {

	struct ModelReduced : ModelParameterized {
		AffineParameterMap affine;
		ModelRBF model;

		ModelReduced();
		~ModelReduced();
		void Create( uint16_t dim, uint8_t paramDim, uint16_t nRBFs );
		void Destroy();
		ModelRBF ComputeModelRBF( const VectorXd& parameter );
		VectorXd VectorField( const VectorXd& x, const VectorXd& parameter );
		MatrixXd Partials( const VectorXd& x, const VectorXd& parameter );
		void OutputText( const char* filename ) const;
	};
}

#endif




