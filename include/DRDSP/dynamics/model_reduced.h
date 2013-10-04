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
		uint32_t dimension, numRBFs;

		ModelReduced();
		~ModelReduced();
		VectorXd Evaluate( const VectorXd &x, const VectorXd &b ) const;
		
		void LoadCentres(char *name);
		void Output(char *name);
		void Create( uint32_t numRB );
		void Destroy();
	};

}

#endif




