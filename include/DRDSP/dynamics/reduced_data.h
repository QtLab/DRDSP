#ifndef INCLUDED_DYNAMICS_REDUCED_DATA
#define INCLUDED_DYNAMICS_REDUCED_DATA
#include "../types.h"
#include "model_orig.h"
#include "../data/data_set.h"

namespace DRDSP {

	struct ReducedData {
		VectorXd *points, *vectors;
		MatrixXd *derivatives;
		uint32_t dimension, count;

		ReducedData();
		ReducedData( uint32_t dim, uint32_t numPoints );
		ReducedData( const ReducedData& rhs );
		ReducedData( ReducedData&& rhs );
		~ReducedData();
		ReducedData& operator=( const ReducedData& rhs );
		ReducedData& operator=( ReducedData&& rhs );
		void Create( uint32_t dim, uint32_t numPoints );
		void Destroy();
		void ComputeData( ModelOriginal& model, const DataSet& data, const VectorXd& parameter, const MatrixXd& W );
		void ComputeBoundingBox( VectorXd& bMin, VectorXd& bMax ) const;
		void WriteData( const char* filename ) const;
		bool ReadData( const char* filename );
	};

}

#endif
