#ifndef INCLUDED_DYNAMICS_REDUCED_DATA_SYSTEM
#define INCLUDED_DYNAMICS_REDUCED_DATA_SYSTEM
#include "reduced_data.h"
#include "../data/data_system.h"

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

		template<typename Family>
		void ComputeDataEmbedded( Family&& family, const DataSystem& data, const MatrixXd& W ) {
			Create( data.numParameters );
			for(uint32_t i=0;i<numParameters;i++) {
				reducedData[i].ComputeDataEmbedded( family(data.parameters[i]), data.dataSets[i], W );
			}
		}

		template<typename Family>
		void ComputeDataCW( Family&& family, const DataSystem& data, const MatrixXd& W ) {
			Create( data.numParameters );
			for(uint32_t i=0;i<numParameters;i++) {
				reducedData[i].ComputeDataCW( family(data.parameters[i]), data.dataSets[i], W );
			}
		}

		template<typename Family>
		void ComputeDataEmbeddedCW( Family&& family, const DataSystem& data, const MatrixXd& W ) {
			Create( data.numParameters );
			for(uint32_t i=0;i<numParameters;i++) {
				reducedData[i].ComputeDataEmbeddedCW( family(data.parameters[i]), data.dataSets[i], W );
			}
		}

		AABB ComputeBoundingBox() const;
		void WritePointsCSV( const char* filePrefix, const char* fileSuffix ) const;
		void WriteVectorsCSV( const char* filePrefix, const char* fileSuffix ) const;
	};

}

#endif
