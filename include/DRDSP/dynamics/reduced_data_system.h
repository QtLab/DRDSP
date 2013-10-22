#ifndef INCLUDED_DYNAMICS_REDUCED_DATA_SYSTEM
#define INCLUDED_DYNAMICS_REDUCED_DATA_SYSTEM
#include "reduced_data.h"

namespace DRDSP {

	struct ReducedDataSystem {
		ReducedData* reducedData;
		uint16_t numParameters;

		ReducedDataSystem();
		ReducedDataSystem( uint16_t N );
		ReducedDataSystem( const ReducedDataSystem& rhs );
		ReducedDataSystem( ReducedDataSystem&& rhs );
		~ReducedDataSystem();
		ReducedDataSystem& operator=( const ReducedDataSystem& rhs );
		ReducedDataSystem& operator=( ReducedDataSystem&& rhs );
		void Create( uint16_t N );
		void Destroy();
		void ComputeData( ModelOriginal& model, const DataSystem& data, const MatrixXd& W );
		AABB ComputeBoundingBox() const;
	};

}

#endif
