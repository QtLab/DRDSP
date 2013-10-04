#ifndef INCLUDED_GEOMETRY_GRASSMANNIAN
#define INCLUDED_GEOMETRY_GRASSMANNIAN
#include "geodesic.h"
#include <Eigen/Dense>

using namespace Eigen;

namespace DRDSP {

	MatrixXd HorizontalComponent( const MatrixXd& W, const MatrixXd& V ) {
		return V - W * (W.adjoint() * V);
	}

	MatrixXd VerticalComponent( const MatrixXd& W, const MatrixXd& V ) {
		return W * (W.adjoint() * V);
	}

	struct Grassmannian : Geodesic<MatrixXd> {
		void Set( const MatrixXd& point, const MatrixXd& tangent );
		MatrixXd Evaluate( double t );
		MatrixXd ParallelTranslate( const MatrixXd &V, double t );
	protected:
		JacobiSVD<MatrixXd> svd;
		typedef JacobiSVD<MatrixXd>::SingularValuesType svType;
	};

}

#endif

