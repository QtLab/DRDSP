#ifndef INCLUDED_DYNAMICS_MONODROMY
#define INCLUDED_DYNAMICS_MONODROMY
#include "../types.h"
#include <vector>
#include <Eigen/Eigenvalues>
#include <limits>

using namespace std;

namespace DRDSP {

	template<typename T>
	T SecantRoot( T x1, T x2, T y1, T y2 ) {
		return x1 - y1 * ( x2 - x1 ) / ( y2 - y1 );
	}

	template<typename V,typename T>
	T DetectPeriod( const vector<V>& orbit, const V& initialVel, T dt, T tolerance = 0.001 ) {
		const auto& initial = orbit[0];
		V initialDir = initialVel.normalized();
		T currentValue = T(0);
		T lastValue = T(0);
		for(size_t i=1;i<orbit.size();++i) {
			lastValue = currentValue;
			currentValue = ( orbit[i] - initial ).dot( initialDir );
			if( lastValue < T(0) && currentValue > T(0) ) {
				T tFactor = SecantRoot( T(0), T(1), lastValue, currentValue );
				auto x = orbit[i-1] * ( T(1) - tFactor ) + orbit[i] * tFactor;
				if( ( x - initial ).squaredNorm() < tolerance*tolerance ) {
					return ( T(i-1) + tFactor ) * dt;
				}
			}
		}
		return T(0);
	}

	template<typename Model>
	MatrixXd ComputeMonodromy( const Model& model, const vector<VectorXd>& samples, double dt ) {
		MatrixXd M;
		M.setIdentity( model.dimension, model.dimension );

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
