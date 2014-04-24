#ifndef INCLUDED_DATA_SECANTS
#define INCLUDED_DATA_SECANTS
#include "../types.h"
#include "data_set.h"
#include <cmath>

namespace DRDSP {

	//! Integer type of weights/counts produced by secant culling
	typedef uint16_t weightType;

	struct SecantsPreComputed;

	/*!
	 * \brief Base class/interface for a set of secants
	 */
	struct Secants {
		size_t count;       //!< Number of secants in the data set
		uint32_t dimension; //!< Dimension of the space

		Secants();
		virtual VectorXd GetSecant( size_t k ) const = 0;                 //!< Get the kth secant (normalized)
		virtual VectorXd GetSecantNoNormalize( size_t k ) const = 0;      //!< Get the kth secant (not necessarily normalized)
		virtual SecantsPreComputed CullSecants( double tolerance ) const;
		SecantsPreComputed CullSecantsDegrees( double degrees ) const;
		SecantsPreComputed CullSecantsRadians( double radians ) const;
	};

	/*!
	 * \brief Set of secants that have been pre-computed
	 */
	struct SecantsPreComputed : Secants {
		vector<VectorXd> secants;       //!< Array of pre-computed unit secants
		vector<weightType> weights;     //!< Array of secant weights produced by culling
		
		void ComputeFromData( const DataSet& dataSet ); //!< Compute this set of secants from the given data set
		VectorXd GetSecant( size_t k ) const;
		VectorXd GetSecantNoNormalize( size_t k ) const;
		SecantsPreComputed CullSecants( double tolerance ) const;
	};

	/*!
	 * \brief Set of secants that are computed on demand from a given DataSet
	 *
	 * For large DataSets, it's not feasible to store all of its secants.
	 * This class will compute each secant on demand.
	 */
	struct SecantsData : Secants {
		const DataSet* data; //!< The DataSet generating the secants

		void SetData( const DataSet& dataSet );          //!< Associate this object with the given data set
		VectorXd GetSecant( size_t k ) const;            //!< Compute and return the kth secant
		VectorXd GetSecantNoNormalize( size_t k ) const; //!< Compute and return the kth secant without normalizing
		static size_t GetIndexI( size_t k, size_t N );
		static size_t GetIndexJ( size_t k, size_t i, size_t N );
	};

}

#endif

