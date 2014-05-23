#ifndef INCLUDED_DATA_SECANTS
#define INCLUDED_DATA_SECANTS
#include "../types.h"
#include "data_system.h"
#include <future>
#include <algorithm>
#include "../geometry/angle.h"

using namespace std;

namespace DRDSP {

	//! Integer type of weights/counts produced by secant culling
	typedef uint16_t weightType;

	/*!
	 * \brief Set of unit secants that have been pre-computed
	 */
	struct Secants {
		size_t count;                //!< Number of secants in the data set
		uint32_t dimension;          //!< Dimension of the space
		vector<VectorXd> secants;    //!< Array of pre-computed unit secants
		vector<weightType> weights;  //!< Array of secant weights produced by culling

		Secants();
		Secants CullSecants( double tolerance ) const;
		Secants CullSecantsAngle( Degreesd degrees ) const;
		Secants CullSecantsAngle( Radiansd radians ) const;
		Secants& ComputeFromData( const DataSet& dataSet ); //!< Compute this set of secants from the given data set
		Secants& ComputeFromData( const DataSet& dataSet, double tolerance ); //!< Compute this set of secants from the given data set with culling
		VectorXd GetSecant( size_t k ) const;
	};

	vector<Secants> ComputeSecants( const DataSystem& data );

	vector<Secants> ComputeSecants( const DataSystem& data, uint32_t numThreads );

	vector<Secants> CullSecants( const vector<Secants>& secants, Degreesd degrees );

	vector<Secants> CullSecants( const vector<Secants>& secants, Degreesd degrees, uint32_t numThreads );

	vector<Secants> ComputeSecants( const DataSystem& data, Degreesd degrees, uint32_t numThreads );
}

#endif
