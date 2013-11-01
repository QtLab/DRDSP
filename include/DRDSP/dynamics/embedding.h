#ifndef INCLUDED_DYNAMICS_EMBEDDING
#define INCLUDED_DYNAMICS_EMBEDDING
/*
|
| Embedding
|
*/
#include "../types.h"
#include "../data/data_set.h"

namespace DRDSP {

	struct Embedding {
		uint32_t oDim, eDim;

		Embedding( uint32_t origDim, uint32_t embedDim );
		virtual VectorXd Evaluate( const VectorXd &x ) const;
		virtual MatrixXd Derivative( const VectorXd &x ) const;
		virtual MatrixXd DerivativeAdjoint( const VectorXd &x ) const;
		virtual MatrixXd Derivative2( const VectorXd &x, uint32_t i ) const;
		MatrixXd ComputeInducedMetric( const VectorXd &x ) const;
		DataSet EmbedData( const DataSet& data ) const;
		DataSystem EmbedData( const DataSystem& data ) const;
	};

	struct EmbeddingCW {
		uint32_t oDim, eDim;

		EmbeddingCW( uint32_t origDim, uint32_t embedDim );
		virtual VectorXd Evaluate( const VectorXd &x ) const;
		virtual double Derivative( const VectorXd &x, uint32_t i, uint32_t j ) const;
		virtual double DerivativeAdjoint( const VectorXd &x, uint32_t i, uint32_t j ) const;
		virtual double Derivative2( const VectorXd &x, uint32_t i, uint32_t j, uint32_t k ) const;
		double ComputeInducedMetric( const VectorXd &x, uint32_t i, uint32_t j ) const;
		DataSet EmbedData( const DataSet& data ) const;
		DataSystem EmbedData( const DataSystem& data ) const;
	};

}

#endif
