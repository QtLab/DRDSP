#ifndef INCLUDED_DYNAMICS_EMBEDDING
#define INCLUDED_DYNAMICS_EMBEDDING
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

}

#endif

