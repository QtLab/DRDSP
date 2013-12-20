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
		uint32_t count,     //! Number of secants in the data set
			     dimension; //! Dimension of the space

		Secants();
		virtual VectorXd GetSecant( uint32_t k ) const = 0; //! Get the kth secant (normalized)
		virtual VectorXd GetSecantNoNormalize( uint32_t k ) const = 0; //! Get the kth secant (not necessarily normalized)
		virtual SecantsPreComputed CullSecants( double tolerance ) const;
		SecantsPreComputed CullSecantsDegrees( double degrees ) const;
		SecantsPreComputed CullSecantsRadians( double radians ) const;
	};

	/*!
	 * \brief Set of secants that have been pre-computed
	 */
	struct SecantsPreComputed : Secants {
		VectorXd* secants;       //! Array of pre-computed unit secants
		weightType* weights;     //! Array of secant weights produced by culling
		
		SecantsPreComputed();
		SecantsPreComputed( const SecantsPreComputed& rhs );
		SecantsPreComputed( SecantsPreComputed&& rhs );
		~SecantsPreComputed();
		SecantsPreComputed& operator=( const SecantsPreComputed& rhs );
		SecantsPreComputed& operator=( SecantsPreComputed&& rhs );
		void ComputeFromData( const DataSet& dataSet ); //! Compute this set of secants from the given data set
		VectorXd GetSecant( uint32_t k ) const;
		VectorXd GetSecantNoNormalize( uint32_t k ) const;
		SecantsPreComputed CullSecants( double tolerance ) const;
	};

	/*!
	 * \brief Set of secants that are computed on demand from a given DataSet
	 *
	 * For large DataSets, it's not feasible to store all of its secants.
	 * This class will compute each secant on demand.
	 */
	// 
	struct SecantsData : Secants {
		const DataSet* data; //! The DataSet generating the secants

		void SetData( const DataSet& dataSet ); //! Associate this object with the given data set
		VectorXd GetSecant( uint32_t k ) const; //! Compute and return the kth secant
		VectorXd GetSecantNoNormalize( uint32_t k ) const; //! Compute and return the kth secant without normalizing
		static uint32_t GetIndexI( uint32_t k, uint32_t N );
		static uint32_t GetIndexJ( uint32_t k, uint32_t i, uint32_t N );
	};

}

#endif

