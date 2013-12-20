#ifndef INCLUDED_DATA_DATA_SET
#define INCLUDED_DATA_DATA_SET
#include "../types.h"

namespace DRDSP {
	
	//! Whether the data points have a type or are untyped
	enum data_t { UNTYPED, TYPED };

	/*!
	 * \brief A set of data points
	 */
	struct DataSet {
		VectorXd* points;   //! Array of data vectors
		uint32_t dimension, //! Dimension of the space
			     count;     //! Number of points in the set
		
		DataSet();
		DataSet( uint32_t dim );
		DataSet( const DataSet& rhs );
		DataSet( DataSet&& rhs );
		~DataSet();
		DataSet& operator=( const DataSet& rhs );
		DataSet& operator=( DataSet&& rhs );
		void Create( uint32_t numPoints, uint32_t dim );
		void Destroy();
		bool LoadBinary( const char* filename );        //! Load a data set from file in binary format
		bool LoadText( const char* filename );          //! Load a data set from file in text format (space separated)
		void WriteCSV( const char* filename ) const;    //! Write the data set to file in CSV format
		DataSet ProjectData( const MatrixXd& W ) const; //! Apply the given projection to this data set
		
		VectorXd& operator[]( uint32_t i ) {
			return points[i];
		}
		
		const VectorXd& operator[]( uint32_t i ) const {
			return points[i];
		}
	};

	/*!
	 * \brief A family of data sets with corresponding parameter values
	 */
	struct DataSystem {
		DataSet*  dataSets;           //! Array of data sets
		VectorXd* parameters;         //! Array of parameter values
		uint32_t  dimension,          //! Dimension of the space in which the data points live
				  maxPoints;          //! Upper limit on the number of data points to read from file
		uint16_t  numParameters;      //! Number of parameters/data sets in the family
		uint8_t   parameterDimension; //! Dimension of the parameter space
		
		DataSystem();
		DataSystem( const DataSystem& rhs );
		DataSystem( DataSystem&& rhs );
		~DataSystem();
		DataSystem& operator=( const DataSystem& rhs );
		DataSystem& operator=( DataSystem&& rhs );
		void Create( uint32_t dim, uint16_t numParams, uint8_t paramDim );
		void Destroy();
		bool Load( const char* filename );                                             //! Load a data system from a text file
		bool LoadSetBinary( const char* filename, uint32_t i );                        //! Load the ith data set from the given file in binary format
		bool LoadSetText( const char* filename, uint32_t i );                          //! Load the ith data set from the given file in text format (space separated)
		void WriteDataSetsCSV( const char* filePrefix, const char* fileSuffix ) const; //! Write the data sets to numbered files with given prefix and suffix
		DataSystem ProjectData( const MatrixXd& W ) const;                             //! Apply the given projection to this data system
	};

}

#endif
