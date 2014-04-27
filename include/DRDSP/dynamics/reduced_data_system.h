#ifndef INCLUDED_DYNAMICS_REDUCED_DATA_SYSTEM
#define INCLUDED_DYNAMICS_REDUCED_DATA_SYSTEM
#include "reduced_data.h"

namespace DRDSP {

	struct ReducedDataSystem {
		vector<ReducedData> reducedData;
		uint32_t numParameters;

		ReducedDataSystem();
		ReducedDataSystem( uint32_t N );
		void Create( uint32_t N );

		template<typename Family>
		void ComputeData( Family&& family, const DataSystem& data, const MatrixXd& W ) {
			Create( data.numParameters );
			for(uint32_t i=0;i<numParameters;i++) {
				reducedData[i].ComputeData( family(data.parameters[i]), data.dataSets[i], W );
			}
		}

		AABB ComputeBoundingBox() const;
		void WritePointsCSV( const char* filePrefix, const char* fileSuffix ) const;
		void WriteVectorsCSV( const char* filePrefix, const char* fileSuffix ) const;
	};

}

#endif
