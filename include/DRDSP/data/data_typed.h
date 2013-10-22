#ifndef INCLUDED_DATA_TYPED
#define INCLUDED_DATA_TYPED
#include "data_set.h"

namespace DRDSP {
	
	struct DataSetTyped : DataSet {
		uint16_t* types;
		
		DataSetTyped();
		DataSetTyped( const DataSetTyped& rhs );
		DataSetTyped( const DataSetTyped&& rhs );
		~DataSetTyped();
		void Create(uint32_t numPoints, uint32_t dim);
		void Destroy();
	};

	struct DataSystemTyped {
		DataSetTyped* dataSets;
		VectorXd* parameters;
		uint32_t dimension,
				 maxPoints;
		uint16_t numParameters;
		uint8_t parameterDimension;
		bool binary;

		DataSystemTyped();
		DataSystemTyped( const DataSystemTyped& rhs );
		DataSystemTyped( const DataSystemTyped&& rhs );
		~DataSystemTyped();
		void Create( uint32_t dim, uint16_t numParams, uint8_t paramDim );
		void Destroy();
		bool Load( const char* filename );
		bool LoadSetBinary( const char* filename, uint32_t i );
		bool LoadSetText( const char* filename, uint32_t i );
	};

}

#endif
