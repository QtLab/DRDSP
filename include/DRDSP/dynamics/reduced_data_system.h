#ifndef INCLUDED_DYNAMICS_REDUCED_DATA_SYSTEM
#define INCLUDED_DYNAMICS_REDUCED_DATA_SYSTEM
#include "../types.h"
#include "model_orig.h"
#include "../data/data_set.h"
#include "reduced_data.h"

namespace DRDSP {

	struct ReducedDataSystem {
		ReducedData* reducedData;
		uint32_t numParameters;

		ReducedDataSystem();
		ReducedDataSystem( uint32_t N );
		ReducedDataSystem( const ReducedDataSystem& rhs );
		ReducedDataSystem( ReducedDataSystem&& rhs );
		~ReducedDataSystem();
		ReducedDataSystem& operator=( const ReducedDataSystem& rhs );
		ReducedDataSystem& operator=( ReducedDataSystem&& rhs );
		void Create( uint32_t N );
		void Destroy();
		void ComputeData( ModelOriginal& model, const DataSystem& data, const MatrixXd& W );
		void ComputeBoundingBox( VectorXd& bMin, VectorXd& bMax ) const;
	};

}

#endif
