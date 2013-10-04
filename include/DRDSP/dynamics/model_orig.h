#ifndef INCLUDED_MODEL_ORIG
#define INCLUDED_MODEL_ORIG
#include "../types.h"

namespace DRDSP {

	struct Embedding {
		uint32_t oDim, eDim;

		Embedding( uint32_t origDim, uint32_t embedDim ) : oDim(origDim), eDim(embedDim) {}

		virtual VectorXd Evaluate( const VectorXd &x ) const {
			return x;
		}

		virtual VectorXd Inverse( const VectorXd &x ) const {
			return x;
		}

		virtual MatrixXd Derivative( const VectorXd &x ) const {
			return MatrixXd::Identity(eDim,oDim);
		}

		virtual double DerivativeCW( const VectorXd &x, uint32_t i, uint32_t j ) const {
			return (i==j)?1.0:0.0;
		}

		virtual MatrixXd Derivative2( const VectorXd &x, uint32_t mu ) const {
			return MatrixXd::Identity(oDim,oDim);
		}

		virtual MatrixXd InverseDerivative( const VectorXd &x ) const {
			return MatrixXd::Identity(oDim,eDim);
		}

		virtual double InverseDerivativeCW( const VectorXd &x, uint32_t i, uint32_t j ) const {
			return (i==j)?1.0:0.0;
		}
	
		virtual MatrixXd InducedMetric( const VectorXd &x ) const {
			MatrixXd aS = Derivative(x);
			return aS.transpose() * aS;
		}

		virtual double InducedMetricCW( const VectorXd &x, uint32_t i, uint32_t j ) const {
			double r = 0.0;
			if( i==j ) {
				double temp;
				for(uint32_t a=0;a<eDim;a++) {
					temp = DerivativeCW(x,a,i);
					r += temp * temp;
				}
			} else {
				for(uint32_t a=0;a<eDim;a++)
					r += DerivativeCW(x,a,i) * DerivativeCW(x,a,j);
			}
			return r;

		}

	};

	struct ModelOriginal {
		uint32_t	pDim,	// dimension of the parameter space
					oDim,	// dimension of the original state space
					eDim;	// dimension of the emedding space

		bool	Euclidean,			// true if the original state space is Euclidean, hence no embedding necessary
				componentWise;		// true if we should evalute the derivatives component-wise, rather than storing a large matrix.

		Embedding* embedding;

		ModelOriginal( uint32_t origDim, uint32_t paramDim ) : oDim(origDim), eDim(origDim), pDim(paramDim), Euclidean(true), componentWise(false), embedding(nullptr) {}
		virtual ~ModelOriginal(){}
		virtual VectorXd VectorField( const VectorXd &x, const VectorXd &b ) = 0;
		virtual double VectorFieldCW( const VectorXd &x, const VectorXd &b, uint32_t i ) { return 0.0; }

		virtual MatrixXd VectorFieldD( const VectorXd &x, const VectorXd &b ) = 0;
		virtual double VectorFieldDCW( const VectorXd &x, const VectorXd &b, uint32_t i, uint32_t j ) { return 0.0; }

		virtual void PrepareOptimizations( const VectorXd &x, const VectorXd &b ) {}
	
	};

}

#endif

