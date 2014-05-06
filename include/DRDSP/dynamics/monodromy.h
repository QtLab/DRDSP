#ifndef INCLUDED_DYNAMICS_MONODROMY
#define INCLUDED_DYNAMICS_MONODROMY
#include "../types.h"
#include <vector>
#include <Eigen/Eigenvalues>

using namespace std;

namespace DRDSP {

	template<typename Model>
	MatrixXd ComputeMonodromy( Model&& model, const vector<VectorXd>& samples, double dt ) {
		MatrixXd M;
		M.setIdentity( model.dimension, model.dimension );

		for( const auto& x : samples ) {
			MatrixXd temp = model.Partials( x ) * M * dt;
			M += temp;
		}
		return M;
	}

	template<typename Model>
	EigenSolver<MatrixXd>::EigenvalueType ComputeFloquetMultipliers( Model&& model, const vector<VectorXd>& samples, double dt ) {
		return EigenSolver<MatrixXd>(
			ComputeMonodromy(
				std::forward<Model>( model ),
				samples,
				dt
			),
			false
		).eigenvalues();
	}

}

#endif
