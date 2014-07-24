#ifndef INCLUDED_AUTO_DIFF
#define INCLUDED_AUTO_DIFF
#include "eigen_dual.h"
#include <Eigen/SparseCore>

namespace DRDSP {

	template<typename F>
	double AutoDerivative( F&& f, double x ) {
		return dual_part( forward<F>(f)( duald(x,1.0) ) );
	}

	template<typename F>
	MatrixXd AutoDerivative( F&& f, const VectorXd& x ) {
		DualMatrixXd df;
		df.setZero( x.size(), x.size() );

		DualVectorXd r = x.cast<duald>();

		for(int64_t j=0;j<df.cols();++j) {
			r[j].y = 1.0;
			df.col( j ) = forward<F>(f)( r );
			r[j].y = 0.0;
		}
		return DualPart( df );
	}

	template<typename F>
	SparseMatrix<double> AutoDerivativeSparse( F&& f, const VectorXd& x ) {
		DualVectorXd r = x.cast<duald>();
		vector<Triplet<double>> triplets;
		triplets.reserve( x.size() );
		VectorXd temp( x.size() );

		for(int64_t j=0;j<x.size();++j) {
			r[j].y = 1.0;
			temp = DualPart( forward<F>(f)( r ) );
			for(int64_t i=0;i<temp.size();++i) {
				if( temp[i] == 0.0 ) continue;
				triplets.emplace_back( i, j, temp[i] );
			}
			r[j].y = 0.0;
		}
		SparseMatrix<double> df((int)x.size(),(int)x.size());
		df.setFromTriplets( cbegin(triplets), cend(triplets) );
		return df;
	}

}

#endif
