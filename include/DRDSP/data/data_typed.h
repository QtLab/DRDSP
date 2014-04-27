#ifndef INCLUDED_DATA_TYPED
#define INCLUDED_DATA_TYPED
#include "data_set.h"

namespace DRDSP {

	//! Whether the data points have a type or are untyped
	enum data_t { UNTYPED, TYPED };

	 /*!
	 * \brief A data set with a corrsponding set of types
	 */
	struct DataSetTyped : DataSet {
		vector<uint16_t> types;
		
		DataSetTyped() = default;
		DataSetTyped( uint32_t numPoints, uint32_t dim );
	};

	/*!
	 * \brief A family of typed data sets with a corresponding set of parameter values
	 */
	struct DataSystemTyped {
		vector<DataSetTyped> dataSets; //!< Array of typed data sets
		vector<VectorXd> parameters;   //!< Array of parameters
		uint32_t dimension;            //!< Dimension of the space in which the data points live
		uint16_t numParameters;        //!< Number of parameters/data sets in the family
		uint8_t parameterDimension;    //!< Dimension of the parameter space


		DataSystemTyped();
		DataSystemTyped( uint32_t dim, uint16_t numParams, uint8_t paramDim );
		bool Load( const char* filename, uint32_t maxPoints = 0 ); //!< Read a data system described by the given text file
		bool LoadSetBinary( const char* filename, uint32_t i, uint32_t maxPoints = 0 ); //!< Load the ith data set from the given file in binary format
		bool LoadSetText( const char* filename, uint32_t i, uint32_t maxPoints = 0 ); //!< Load the ith data set from the given file in text format (space separated)
	};

}

#endif
