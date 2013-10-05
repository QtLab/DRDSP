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
		RadialFunction* rbfs;
		uint32_t dimension, parameterDimension, numRBFs;

		ModelReduced();
		~ModelReduced();
		void Create( uint32_t dim, uint32_t paramDim, uint32_t nRBFs );
		void Destroy();
		VectorXd Evaluate( const VectorXd &x, const VectorXd &parameter ) const;
		void OutputText( const char *filename ) const;
	};
}

#endif




