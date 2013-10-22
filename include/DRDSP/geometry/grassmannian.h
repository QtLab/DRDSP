#ifndef INCLUDED_GEOMETRY_GRASSMANNIAN
#define INCLUDED_GEOMETRY_GRASSMANNIAN
#include "geodesic.h"
#include <Eigen/Dense>

using namespace Eigen;

namespace DRDSP {

	struct Grassmannian : Geodesic<MatrixXd> {
		void Set( const MatrixXd& point, const MatrixXd& tangent );
		MatrixXd Evaluate( double t );
		MatrixXd ParallelTranslate( const MatrixXd &V, double t );
		static MatrixXd HorizontalComponent( const MatrixXd& W, const MatrixXd& V );
		static MatrixXd VerticalComponent( const MatrixXd& W, const MatrixXd& V );
	protected:
		JacobiSVD<MatrixXd> svd;
		typedef JacobiSVD<MatrixXd>::SingularValuesType svType;
	};

}

#endif

