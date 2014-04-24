#ifndef INCLUDED_DATA_DATA_SET
#define INCLUDED_DATA_DATA_SET
#include "../types.h"
#include <vector>
#include <cmath>

using namespace std;

namespace DRDSP {
	
	//! Whether the data points have a type or are untyped
	enum data_t { UNTYPED, TYPED };

	/*!
	 * \brief A set of data points
	 */
	struct DataSet {
		vector<VectorXd> points;   //!< Array of data vectors
		uint32_t dimension; //!< Dimension of the space
		
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

	/*!
	 * \brief A family of data sets with corresponding parameter values
	 */
	struct DataSystem {
		vector<DataSet> dataSets;     //!< Array of data sets
		vector<VectorXd> parameters;  //!< Array of parameter values
		uint32_t dimension,           //!< Dimension of the space in which the data points live
				 maxPoints;           //!< Upper limit on the number of data points to read from file
		uint16_t numParameters;       //!< Number of parameters/data sets in the family
		uint8_t parameterDimension;   //!< Dimension of the parameter space
		
		DataSystem();
		DataSystem( uint32_t dim, uint16_t numParams, uint8_t paramDim );
		bool Load( const char* filename );                                             //!< Load a data system from a text file
		bool LoadSetBinary( const char* filename, uint32_t i );                        //!< Load the ith data set from the given file in binary format
		bool LoadSetText( const char* filename, uint32_t i );                          //!< Load the ith data set from the given file in text format (space separated)
		void WriteDataSetsCSV( const char* filePrefix, const char* fileSuffix ) const; //!< Write the data sets to numbered files with given prefix and suffix
		DataSystem ProjectData( const MatrixXd& W ) const;                             //!< Apply the given projection to this data system
	};

	vector<double> NormDifferences( const vector<VectorXd>& s1, const vector<VectorXd>& s2 );
	double RMSDifference( const vector<VectorXd>& s1, const vector<VectorXd>& s2 );

}

#endif
