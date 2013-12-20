#ifndef INCLUDED_DYNAMICS_EMBEDDING
#define INCLUDED_DYNAMICS_EMBEDDING
#include "../types.h"
#include "../data/data_set.h"

namespace DRDSP {

	/*!
	 * \brief An embedding of a state space into R^n
	 */
	struct Embedding {
		uint32_t oDim, //! Dimension of the original state space
			     eDim; //! Dimension of the embedding space, R^n

		Embedding( uint32_t origDim, uint32_t embedDim );
		virtual VectorXd Evaluate( const VectorXd &x ) const; //! Embeds the given state
		virtual MatrixXd Derivative( const VectorXd &x ) const; //! The partial derivatives of the embedding with respect to the state
		virtual MatrixXd DerivativeAdjoint( const VectorXd &x ) const; //! The linear algebraic adjoint of the derivative
		virtual MatrixXd Derivative2( const VectorXd &x, uint32_t i ) const; //! The second partial derivatives of the ith component with respect to the state
		MatrixXd ComputeInducedMetric( const VectorXd &x ) const; //! The pullback metric on the state space from R^n
		DataSet EmbedData( const DataSet& data ) const; //! Apply this embedding to the data set
		DataSystem EmbedData( const DataSystem& data ) const; //! Apply this embedding to the data system
	};

	/*!
	 * \brief An embedding of a state space into R^n
	 * whose matrices are to be evaluated component-wise
	 *
	 * This version of the embedding is for very high-dimensional systems
	 * whose matrix of partial derivatives is too large to store in memory.
	 * These systems will have their derivatives evaluated one element at a time.
	 */
	struct EmbeddingCW {
		uint32_t oDim, //! Dimension of the original state space
			     eDim; //! Dimension of the embedding space, R^n

		EmbeddingCW( uint32_t origDim, uint32_t embedDim );
		virtual VectorXd Evaluate( const VectorXd &x ) const; //! Embeds the given state
		virtual double Derivative( const VectorXd &x, uint32_t i, uint32_t j ) const; //! The partial derivatives of the ith component of the embedding with respect to the jth component of the state
		virtual double DerivativeAdjoint( const VectorXd &x, uint32_t i, uint32_t j ) const; //! The (i,j) component of the linear algebraic adjoint of the derivative
		virtual double Derivative2( const VectorXd &x, uint32_t i, uint32_t j, uint32_t k ) const; //! The second partial derivatives of the ith component with respect to the j and k components of the state
		double ComputeInducedMetric( const VectorXd &x, uint32_t i, uint32_t j ) const; //! The (i,j) component of the pullback metric on the state space from R^n
		DataSet EmbedData( const DataSet& data ) const; //! Apply this embedding to the data set
		DataSystem EmbedData( const DataSystem& data ) const; //! Apply this embedding to the data system
	};

}

#endif
