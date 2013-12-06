#ifndef INCLUDED_DATA_SECANTS
#define INCLUDED_DATA_SECANTS
#include "../types.h"
#include "data_set.h"
#include <cmath>

namespace DRDSP {

	typedef uint16_t weightType;

	struct SecantsPreComputed;

	struct Secants {
		uint32_t count, dimension;

		Secants();
		virtual VectorXd GetSecant( uint32_t k ) const = 0;
		virtual VectorXd GetSecantNoNormalize( uint32_t k ) const = 0;
		virtual SecantsPreComputed CullSecants( double tolerance ) const;
		SecantsPreComputed CullSecantsDegrees( double degrees ) const;
		SecantsPreComputed CullSecantsRadians( double radians ) const;
	};

	struct SecantsPreComputed : Secants {
		VectorXd* secants;       // Pre-computed secants
		weightType* weights;     // Secant weights produced by culling

		SecantsPreComputed();
		SecantsPreComputed( const SecantsPreComputed& rhs );
		SecantsPreComputed( SecantsPreComputed&& rhs );
		~SecantsPreComputed();
		SecantsPreComputed& operator=( const SecantsPreComputed& rhs );
		SecantsPreComputed& operator=( SecantsPreComputed&& rhs );
		void ComputeFromData( const DataSet& dataSet );
		VectorXd GetSecant( uint32_t k ) const;
		VectorXd GetSecantNoNormalize( uint32_t k ) const;
		SecantsPreComputed CullSecants( double tolerance ) const;
	};

	// For large sets, we determine each secant from the data points on demand
	struct SecantsData : Secants {
		const DataSet* data;

		void SetData( const DataSet& dataSet );
		VectorXd GetSecant( uint32_t k ) const;
		VectorXd GetSecantNoNormalize( uint32_t k ) const;
		static uint32_t GetIndexI( uint32_t k, uint32_t N );
		static uint32_t GetIndexJ( uint32_t k, uint32_t i, uint32_t N );
	};

}

#endif

