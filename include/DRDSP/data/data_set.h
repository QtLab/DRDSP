#ifndef INCLUDED_DATA_DATA_SET
#define INCLUDED_DATA_DATA_SET
#include "../types.h"
#include <vector>

using namespace std;

namespace DRDSP {

	/*!
	 * \brief A set of data points
	 */
	struct DataSet {
		vector<VectorXd> points;   //!< Array of data vectors
		uint32_t dimension;        //!< Dimension of the space
		
		DataSet() : dimension(0) {}
		DataSet( size_t numPoints, uint32_t dim );
		bool LoadBinary( const char* filename );        //!< Load a data set from file in binary format
		bool LoadText( const char* filename );          //!< Load a data set from file in text format (space separated)
		void WriteCSV( const char* filename ) const;    //!< Write the data set to file in CSV format
		DataSet ProjectData( const MatrixXd& W ) const; //!< Apply the given projection to this data set
		
		VectorXd& operator[]( size_t i ) {
			return points[i];
		}
		
		const VectorXd& operator[]( size_t i ) const {
			return points[i];
		}
	};

	template<typename E>
	DataSet EmbedData( const E& embedding, const DataSet& data ) {
		DataSet r( data.points.size(), embedding.eDim );
		for(uint32_t i=0;i<r.points.size();++i) {
			r[i] = embedding(data[i]);
		}
		return r;
	}

	struct DataComparisonResult {
		double maxDifference, rmsDifference, maxMinDifference;
		vector<double> differences;

		DataComparisonResult( size_t n ) : differences(n) {}
	};

	DataComparisonResult CompareData( const vector<VectorXd>& s1, const vector<VectorXd>& s2 );

}

#endif
