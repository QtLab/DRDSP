#ifndef INCLUDED_MODEL_REDUCED
#define INCLUDED_MODEL_REDUCED
#include <fstream>
#include "../types.h"
#include "affineParameterMap.h"
#include "model_rbf.h"
#include "reduced_data_system.h"

using namespace std;

namespace DRDSP {

	struct ModelReduced {
		AffineParameterMap affine;
		uint16_t dimension;
		uint8_t parameterDimension;
		ModelRBF model;

		ModelReduced();
		~ModelReduced();
		void Create( uint16_t dim, uint8_t paramDim, uint16_t nRBFs );
		void Destroy();
		ModelRBF ComputeModelRBF( const VectorXd& parameter );
		VectorXd Evaluate( const VectorXd& x, const VectorXd& parameter );
		void OutputText( const char* filename ) const;
	};
}

#endif




