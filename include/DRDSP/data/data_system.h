#ifndef INCLUDED_DATA_DATA_SYSTEM
#define INCLUDED_DATA_DATA_SYSTEM
#include "data_set.h"

namespace DRDSP {
	
	/*!
	 * \brief A family of data sets with corresponding parameter values
	 */
	struct DataSystem {
		vector<DataSet> dataSets;     //!< Array of data sets
		vector<VectorXd> parameters; //!< Array of parameter values
		uint32_t dimension,           //!< Dimension of the space in which the data points live
		         numParameters,       //!< Number of parameters/data sets in the family
		         parameterDimension;  //!< Dimension of the parameter space

		DataSystem();

		DataSystem( uint32_t dim, uint32_t numParams, uint32_t paramDim );

		//! Load a data system from a text file
		bool Load( const char* filename, uint32_t maxPoints = 0 );

		//! Load the ith data set from the given file in binary format
		bool LoadSetBinary( const char* filename, uint32_t i, uint32_t maxPoints = 0 );

		//! Load the ith data set from the given file in text format (space separated)
		bool LoadSetText( const char* filename, uint32_t i, uint32_t maxPoints = 0 );

		//! Write the data sets to numbered files with given prefix and suffix
		void WriteDataSetsCSV( const char* filePrefix, const char* fileSuffix ) const;

		//! Apply the given projection to this data system
		DataSystem ProjectData( const MatrixXd& W ) const;
	
	};

	template<typename E>
	DataSystem EmbedData( const E& embedding, const DataSystem& data ) {
		DataSystem r( embedding.eDim, data.numParameters, data.parameterDimension );
		r.parameters = data.parameters;
		
		for(uint32_t i=0;i<r.numParameters;++i) {
			r.dataSets[i].points.resize( data.dataSets[i].points.size() );

			for(uint32_t j=0;j<r.dataSets[i].points.size();++j)
				r.dataSets[i][j] = embedding(data.dataSets[i][j]);
		}

		return r;
	}

}

#endif
