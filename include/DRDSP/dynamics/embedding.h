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

	struct IdentityEmbeddingCW : Embedding {

		IdentityEmbeddingCW( uint32_t dimension ) : Embedding(dimension,dimension) {}

		VectorXd Evaluate( const VectorXd& x ) const {
			return x;
		}

		double Derivative( const VectorXd&, uint32_t i, uint32_t j ) const {
			return (i==j)?1.0:0.0;
		}

		double DerivativeAdjoint( const VectorXd&, uint32_t i, uint32_t j ) const {
			return (i==j)?1.0:0.0;
		}

		double Derivative2( const VectorXd&, uint32_t, uint32_t j, uint32_t k ) const {
			return (j==k)?1.0:0.0;
		}

	};

	template<typename E>
	MatrixXd ComputeInducedMetric( const E& embedding, const VectorXd& x ) {
		MatrixXd deriv = embedding.Derivative(x);
		return deriv.transpose() * deriv;
	}

	template<typename E>
	double ComputeInducedMetric( const E& embedding, const VectorXd& x, uint32_t i, uint32_t j ) {
		double r = 0.0;
		if( i==j ) {
			double temp;
			for(uint32_t a=0;a<eDim;++a) {
				temp = embedding.Derivative(x,a,i);
				r += temp * temp;
			}
		} else {
			for(uint32_t a=0;a<eDim;++a)
				r += embedding.Derivative(x,a,i) * embedding.Derivative(x,a,j);
		}
		return r;
	}

}

#endif
