#ifndef INCLUDED_PROJ_INVERSE
#define INCLUDED_PROJ_INVERSE
#include "../data/data_system.h"
#include "../eigen_affine.h"

namespace DRDSP {

	struct ProjInverse {
		AffineXd vertical;
		MatrixXd W;

		VectorXd operator()( const VectorXd& state, const VectorXd& parameter ) const {
			return W * state + vertical( parameter );
		}

		MatrixXd Derivative( const VectorXd&, const VectorXd& ) const {
			return W;
		}

		void Compute( const MatrixXd& W, const DataSystem& data );

	};

	ProjInverse ComputeInverse( const MatrixXd& W, const DataSystem& data );

	double ComputeInverseCost( const ProjInverse& inverse, const DataSystem& data );

}

#endif
