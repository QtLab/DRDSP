#ifndef INCLUDED_DATA_DATA_SET
#define INCLUDED_DATA_DATA_SET
#include "../types.h"
#include "../geometry/ball.h"
#include <vector>

using namespace std;

namespace DRDSP {

	/**
	 * \brief A set of data points
	 */
	struct DataSet {
		uint32_t dimension = 0;    ///< Dimension of the space
		vector<VectorXd> points;   ///< Array of data vectors
		
		DataSet() = default;
		explicit DataSet( size_t numPoints ) : points(numPoints) {}
		DataSet( size_t numPoints, uint32_t dim );
		bool LoadBinary( const char* filename );         ///< Load a data set from file in binary format
		bool LoadText( const char* filename );           ///< Load a data set from file in text format (space separated)
		void WriteCSV( const char* filename ) const;     ///< Write the data set to file in CSV format
		void WriteBinary( const char* filename ) const;  ///< Write the data set to file in binary format
		
		VectorXd& operator[]( size_t i ) {
			return points[i];
		}
		
		const VectorXd& operator[]( size_t i ) const {
			return points[i];
		}
	};

	template<typename F>
	DataSet Map( const F& f, const DataSet& data ) {
		DataSet r( data.points.size() );
		for(size_t i=0;i<r.points.size();++i) {
			r[i] = f(data[i]);
		}
		r.dimension = (uint32_t)r[0].size();
		return r;
	}

	template<typename Derived>
	DataSet MapMatrix( const MatrixBase<Derived>& A, const DataSet& data ) {
		DataSet r( data.points.size() );
		for(size_t i=0;i<r.points.size();++i) {
			r[i].noalias() = A * data[i];
		}
		r.dimension = (uint32_t)r[0].size();
		return r;
	}

	struct DataComparisonResult {
		double maxDifference, rmsDifference, maxMinDifference;
		vector<double> differences;

		DataComparisonResult( size_t n ) : differences(n) {}
	};

	DataComparisonResult CompareData( const vector<VectorXd>& s1, const vector<VectorXd>& s2 );

	BallXd ComputeBall( const vector<VectorXd>& points );


}

#endif
