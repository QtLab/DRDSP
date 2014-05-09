#ifndef INCLUDED_AUTO_DIFF
#define INCLUDED_AUTO_DIFF
#include "eigen_dual.h"

namespace DRDSP {

	template<typename F>
	double Derivative( F&& f, double x ) {
		return dual_part( f( dual<T>(x,1.0) ) );
	}

	template<typename F>
	MatrixXd Derivative( F&& f, const VectorXd& x ) {
		DualMatrixXd df;
		df.setZero( x.size(), x.size() );

		DualVectorXd r = x.cast<duald>();
		DualVectorXd rj;

		for(int64_t j=0;j<df.cols();++j) {
			rj = r;
			rj[j].y = 1.0;
			df.col( j ) = f( rj );
		}
		return DualPart( df );
	}


}

#endif
