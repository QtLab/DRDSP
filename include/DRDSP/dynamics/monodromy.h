#ifndef INCLUDED_DYNAMICS_MONODROMY
#define INCLUDED_DYNAMICS_MONODROMY
#include "../types.h"
#include <vector>
#include <Eigen/Eigenvalues>

using namespace std;

namespace DRDSP {

	template<typename T>
	T SecantRoot( T x1, T x2, T y1, T y2 ) {
		return x1 - y1 * ( x2 - x1 ) / ( y2 - y1 );
	}

	template<typename V>
	double DetectPeriod( const vector<V>& orbit, const V& initialVel, double dt, double tolerance = 0.001 ) {
		const auto& initial = orbit[0];
		V initialDir = initialVel.normalized();
		double currentValue = 0.0;
		double lastValue = 0.0;
		for(size_t i=1;i<orbit.size();++i) {
			lastValue = currentValue;
			currentValue = ( orbit[i] - initial ).dot( initialDir );
			if( lastValue < 0.0 && currentValue > 0.0 ) {
				double tFactor = SecantRoot( 0.0, 1.0, lastValue, currentValue );
				auto x = orbit[i-1] * ( 1.0 - tFactor ) + orbit[i] * tFactor;
				if( ( x - initial ).squaredNorm() < tolerance*tolerance ) {
					return ( double(i-1) + tFactor ) * dt;
				}
			}
		}
		return 0.0;
	}

	template<typename Model>
	MatrixXd ComputeMonodromy( const Model& model, const vector<VectorXd>& samples, double dt ) {
		MatrixXd M;
		M.setIdentity( model.stateDim, model.stateDim );

		for( const auto& x : samples ) {
			MatrixXd temp = model.Partials( x ) * M * dt;
			M += temp;
		}
		return M;
	}

	template<typename Model>
	EigenSolver<MatrixXd>::EigenvalueType ComputeFloquetMultipliers( const Model& model, const vector<VectorXd>& samples, double dt ) {
		return EigenSolver<MatrixXd>(
			ComputeMonodromy( model, samples, dt ),
			false
		).eigenvalues();
	}

}

#endif
