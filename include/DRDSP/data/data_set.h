#ifndef INCLUDED_DATA_DATA_SET
#define INCLUDED_DATA_DATA_SET
#include "../types.h"

namespace DRDSP {
	
	enum data_t { UNTYPED, TYPED };

	struct DataSet {
		VectorXd* points;
		uint32_t dimension, count;
		
		DataSet();
		DataSet( uint32_t dim );
		DataSet( const DataSet& rhs );
		DataSet( DataSet&& rhs );
		~DataSet();
		DataSet& operator=( const DataSet& rhs );
		DataSet& operator=( DataSet&& rhs );
		void Create( uint32_t numPoints, uint32_t dim );
		void Destroy();
		bool LoadBinary( const char* filename );
		bool LoadText( const char* filename );
		void WriteCSV( const char* filename ) const;
		DataSet ProjectData( const MatrixXd& W ) const;
		
		VectorXd& operator[]( uint32_t i ) {
			return points[i];
		}
		
		const VectorXd& operator[]( uint32_t i ) const {
			return points[i];
		}
	};

	struct DataSystem {
		DataSet*  dataSets;
		VectorXd* parameters;
		uint32_t  dimension,
				  maxPoints;
		uint16_t  numParameters;
		uint8_t   parameterDimension;
		
		DataSystem();
		DataSystem( const DataSystem& rhs );
		DataSystem( DataSystem&& rhs );
		~DataSystem();
		DataSystem& operator=( const DataSystem& rhs );
		DataSystem& operator=( DataSystem&& rhs );
		void Create( uint32_t dim, uint16_t numParams, uint8_t paramDim );
		void Destroy();
		bool Load( const char* filename );
		bool LoadSetBinary( const char* filename, uint32_t i );
		bool LoadSetText( const char* filename, uint32_t i );
		void WriteDataSetsCSV( const char* filePrefix, const char* fileSuffix ) const;
		DataSystem ProjectData( const MatrixXd& W ) const;
	};

}

#endif
