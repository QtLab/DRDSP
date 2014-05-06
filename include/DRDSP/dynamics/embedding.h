#ifndef INCLUDED_DYNAMICS_EMBEDDING
#define INCLUDED_DYNAMICS_EMBEDDING
#include "../types.h"

namespace DRDSP {

	/*!
	 * \brief Base for an embedding of a state space into R^n
	 */
	struct Embedding {
		uint32_t oDim, //!< Dimension of the original state space
			     eDim; //!< Dimension of the embedding space, R^n

		Embedding( uint32_t originalDimension, uint32_t embeddedDimension ) : oDim(originalDimension), eDim(embeddedDimension) {}

	};

	struct IdentityEmbedding : Embedding {

		IdentityEmbedding( uint32_t dimension ) : Embedding(dimension,dimension) {}

		VectorXd operator()( const VectorXd& x ) const {
			return x;
		}

		MatrixXd Derivative( const VectorXd& ) const {
			return MatrixXd::Identity(eDim,oDim);
		}

		MatrixXd DerivativeAdjoint( const VectorXd& ) const {
			return MatrixXd::Identity(oDim,eDim);
		}

		MatrixXd Derivative2( const VectorXd&, uint32_t ) const {
			return MatrixXd::Identity(oDim,oDim);
		}

	};

	template<typename E>
	MatrixXd ComputeInducedMetric( const E& embedding, const VectorXd& x ) {
		MatrixXd deriv = embedding.Derivative(x);
		return deriv.transpose() * deriv;
	}

}

#endif
