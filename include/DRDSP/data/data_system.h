#ifndef INCLUDED_DATA_DATA_SYSTEM
#define INCLUDED_DATA_DATA_SYSTEM
#include "data_set.h"

namespace DRDSP {
	
	/*!
	 * \brief A family of data sets with corresponding parameter values
	 */
	struct DataSystem {
		uint32_t dimension,           //!< Dimension of the space in which the data points live
		         numParameters,       //!< Number of parameters/data sets in the family
		         parameterDimension;  //!< Dimension of the parameter space
		vector<DataSet> dataSets;     //!< Array of data sets
		vector<VectorXd> parameters;  //!< Array of parameter values

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

			for(size_t j=0;j<r.dataSets[i].points.size();++j)
				r.dataSets[i][j] = embedding(data.dataSets[i][j]);
		}

		return r;
	}

	template<typename E>
	DataSystem EmbedData( const E& embedding, const DataSystem& data, uint32_t numThreads ) {
		DataSystem r( embedding.eDim, data.numParameters, data.parameterDimension );
		r.parameters = data.parameters;
		vector<future<void>> futures(numThreads);
		for(uint32_t i=0;i<data.numParameters;i+=numThreads) {
			uint32_t N = min( data.numParameters - i, numThreads );
			for(uint32_t j=0;j<N;++j) {
				futures[j] = async( launch::async,
					[&]( DataSet& embedded, const DataSet& original ) {
						embedded.points.resize( original.points.size() );
						for(size_t k=0;k<embedded.points.size();++k)
							embedded[k] = embedding( original[k] );
					},
					ref(r.dataSets[i+j]), cref(data.dataSets[i+j])
				);
			}
			for(uint32_t j=0;j<N;++j) {
				futures[j].wait();
			}
		}
		return r;
	}

}

#endif
