#ifndef INCLUDED_DATA_DATA_SET
#define INCLUDED_DATA_DATA_SET
#include "../types.h"

namespace DRDSP {
	
	enum data_t { UNTYPED, TYPED };

	struct DataSet {
		VectorXd* points;
		uint32_t dimension, count;
		
		DataSet();
		DataSet( const DataSet& rhs );
		DataSet( DataSet&& rhs );
		~DataSet();
		void Create( uint32_t numPoints, uint32_t dim );
		void Destroy();
		DataSet ProjectData( const MatrixXd& W ) const;
		
		VectorXd& operator[]( uint32_t i ) {
			return points[i];
		}
		
		const VectorXd& operator[]( uint32_t i ) const {
			return points[i];
		}
	};

	struct DataSystem {
		DataSet *dataSets;
		VectorXd *parameters;
		uint32_t numParameters,
		         dimension,
				 parameterDimension,
				 maxPoints;
		
		DataSystem();
		DataSystem( const DataSystem& rhs );
		DataSystem( DataSystem&& rhs );
		~DataSystem();
		void Create( uint32_t numParams, uint32_t paramDim, uint32_t dim );
		void Destroy();
		bool Load( const char* filename );
		bool LoadSetBinary( const char* filename, uint32_t i );
		bool LoadSetText( const char* filename, uint32_t i );
		DataSystem ProjectData( const MatrixXd& W ) const;
	};

}

#endif
