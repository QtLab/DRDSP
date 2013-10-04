#ifndef INCLUDED_DATA_SECANTS
#define INCLUDED_DATA_SECANTS
#include "../types.h"
#include "data_set.h"
#include <cmath>

namespace DRDSP {

	typedef uint16_t weightType;

	struct Secants {
		VectorXd* secants;     // Pre-computed secants
		weightType* weights;     // Secant weights produced by culling
		uint32_t count;        // Number of secants
		uint32_t dimension;
		bool preComputed;
		bool weighted;

		Secants();
		Secants( const Secants& rhs );
		Secants( Secants&& rhs );
		~Secants();
		void ComputeFromData( const DataSet& dataSet, size_t preComputeSize = 0 );
		VectorXd GetSecant( uint32_t k ) const;
	
	protected:
		void PreCompute();

		// For large sets, we determine each secant from the data points on demand
		const DataSet* data;
		static uint32_t GetIndexI( uint32_t k, uint32_t N );
		static uint32_t GetIndexJ( uint32_t k, uint32_t i, uint32_t N );
	};

	struct SecantCulling {
		double tolerance;

		void SetAngularTolerance( double angle ) {
			tolerance = cos( angle );
		}
	
		void SetAngularToleranceDegrees( double angle ) {
			tolerance = cos( angle * (M_PI/180.0) );
		}

		Secants CullSecants( const Secants& input )  const;
	};

}

#endif

