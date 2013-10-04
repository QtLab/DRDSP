#ifndef INCLUDED_DYNAMICS_AFFINEPARAMETERMAP
#define INCLUDED_DYNAMICS_AFFINEPARAMETERMAP
#include "../types.h"

namespace DRDSP {

	struct AffineParameterMap {
		uint32_t M, R, N, dim;
		MatrixXd coeffs;
	
		void Init( uint32_t x, uint32_t y, uint32_t d ) {
			M = x;
			N = y;
			dim = d;
			R = N * (dim+1);
			coeffs.setZero(M,R);
		}

		MatrixXd Evaluate( const VectorXd &x ) const {
			return coeffs * GetLambda(x);
		}
		MatrixXd GetLambda( const VectorXd &x ) const {
			MatrixXd Lambda;
			Lambda.setZero(R,N);
			Lambda.block(0,0,N,N).setIdentity();
			for(uint32_t i=0;i<dim;i++) {
				Lambda.block(N*(i+1),0,N,N).setIdentity();
				Lambda.block(N*(i+1),0,N,N) *= x(i);
			}
			return Lambda;
		}

	};

}

#endif

