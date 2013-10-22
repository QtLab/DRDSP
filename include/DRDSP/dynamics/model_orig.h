#ifndef INCLUDED_DYNAMICS_MODEL_ORIG
#define INCLUDED_DYNAMICS_MODEL_ORIG
#include "embedding.h"

namespace DRDSP {

	struct ModelOriginal {
		uint32_t oDim,  // dimension of the original state space
		         eDim;  // dimension of the emedding space
		uint8_t  pDim;   // dimension of the parameter space

		bool Euclidean,         // true if the original state space is Euclidean, hence no embedding necessary
		     componentWise;     // true if we should evalute the derivatives component-wise, rather than storing a large matrix.

		Embedding* embedding;

		ModelOriginal( uint32_t origDim, uint32_t paramDim ) : oDim(origDim), eDim(origDim), pDim(paramDim), Euclidean(true), componentWise(false), embedding(nullptr) {}
		virtual ~ModelOriginal() {}
		virtual VectorXd VectorField( const VectorXd &x, const VectorXd &b ) = 0;
		virtual double VectorFieldCW( const VectorXd &x, const VectorXd &b, uint32_t i ) { return 0.0; }

		virtual MatrixXd VectorFieldD( const VectorXd &x, const VectorXd &b ) = 0;
		virtual double VectorFieldDCW( const VectorXd &x, const VectorXd &b, uint32_t i, uint32_t j ) { return 0.0; }

		virtual void PrepareOptimizations( const VectorXd &x, const VectorXd &b ) {}
	
	};

}

#endif

