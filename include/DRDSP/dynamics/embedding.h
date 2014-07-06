#ifndef INCLUDED_DYNAMICS_EMBEDDING
#define INCLUDED_DYNAMICS_EMBEDDING
#include "../types.h"

namespace DRDSP {

	/**
	 * \brief Base for an embedding of a state space into R^n
	 */
	struct Embedding {
		uint32_t sourceDim,   ///< Dimension of the original state space
			     embedDim;    ///< Dimension of the embedding space, R^n

		Embedding( uint32_t sourceDim, uint32_t embedDim ) : sourceDim(sourceDim), embedDim(embedDim) {}

	};

	struct IdentityEmbedding : Embedding {

		IdentityEmbedding( uint32_t dim ) : Embedding(dim,dim) {}

		VectorXd operator()( const VectorXd& x ) const {
			return x;
		}

		MatrixXd Derivative( const VectorXd& ) const {
			return MatrixXd::Identity(embedDim,sourceDim);
		}

		MatrixXd DerivativeAdjoint( const VectorXd& ) const {
			return MatrixXd::Identity(sourceDim,embedDim);
		}

		MatrixXd Derivative2( const VectorXd&, uint32_t ) const {
			return MatrixXd::Identity(sourceDim,sourceDim);
		}

	};

	template<typename E>
	MatrixXd ComputeInducedMetric( const E& embedding, const VectorXd& x ) {
		MatrixXd deriv = embedding.Derivative(x);
		return deriv.transpose() * deriv;
	}

}

#endif
