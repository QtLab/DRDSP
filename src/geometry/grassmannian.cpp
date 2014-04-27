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

	svType DH = svd.singularValues() * t;

	MatrixXd Wt = position;
	Wt *= svd.matrixV() * DH.array().cos().matrix().asDiagonal(); 
	Wt += svd.matrixU() * DH.array().sin().matrix().asDiagonal();
	Wt *= svd.matrixV().transpose();

	return Wt;
}

MatrixXd Geodesic::ParallelTranslate( const MatrixXd& V, double t ) const {
	if( t == 0.0 ) return velocity;
	
	svType DH = svd.singularValues() * t;
	MatrixXd UtV = svd.matrixU().transpose() * velocity;

	MatrixXd Vt = position;
	Vt *= -svd.matrixV() * DH.array().sin().matrix().asDiagonal(); 
	Vt += svd.matrixU() * DH.array().cos().matrix().asDiagonal();
	Vt *= UtV;
	Vt += V - svd.matrixU() * UtV;

	return Vt;
}
