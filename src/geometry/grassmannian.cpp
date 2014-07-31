#include <DRDSP/geometry/grassmannian.h>

using namespace DRDSP;
using namespace DRDSP::Grassmannian;

MatrixXd Grassmannian::HorizontalComponent( const MatrixXd& W, const MatrixXd& V ) {
	return V - W * (W.adjoint() * V);
}

MatrixXd Grassmannian::VerticalComponent( const MatrixXd& W, const MatrixXd& V ) {
	return W * (W.adjoint() * V);
}

void Geodesic::Set( const MatrixXd& x, const MatrixXd& v ) {
	position = x;
	velocity = v;
	svd.compute( velocity, ComputeThinU | ComputeThinV );
}

MatrixXd Geodesic::operator()( double t ) const {
	if( t == 0.0 ) return position;

	svType st = svd.singularValues() * t;

	return (
		     position * svd.matrixV() * st.array().cos().matrix().asDiagonal()
	       + svd.matrixU() * st.array().sin().matrix().asDiagonal()
	) * svd.matrixV().transpose();
}

MatrixXd Geodesic::ParallelTranslate( const MatrixXd& V, double t ) const {
	if( t == 0.0 ) return V;
	
	svType st = svd.singularValues() * t;

	return V + ( 
	      svd.matrixU() * st.array().cos().matrix().asDiagonal()
		- position * svd.matrixV() * st.array().sin().matrix().asDiagonal()
	    - svd.matrixU()
	) * svd.matrixU().transpose() * V;
}
