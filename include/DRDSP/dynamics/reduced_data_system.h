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
		void ComputeData( ModelParameterized& model, const DataSystem& data, const MatrixXd& W );
		void ComputeData( ModelParameterizedCW& model, const DataSystem& data, const MatrixXd& W );
		void ComputeData( ModelParameterizedEmbedded& model, const DataSystem& data, const MatrixXd& W );
		void ComputeData( ModelParameterizedEmbeddedCW& model, const DataSystem& data, const MatrixXd& W );
		AABB ComputeBoundingBox() const;
		void WritePointsCSV( const char* filePrefix, const char* fileSuffix ) const;
		void WriteVectorsCSV( const char* filePrefix, const char* fileSuffix ) const;
	};

}

#endif
