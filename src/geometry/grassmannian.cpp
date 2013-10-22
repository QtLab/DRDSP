#include <DRDSP/geometry/grassmannian.h>

using namespace DRDSP;

void Grassmannian::Set( const MatrixXd& point, const MatrixXd& tangent ) {
	Point = point;
	Tangent = tangent;
	svd.compute(Tangent,ComputeThinU | ComputeThinV);
}

MatrixXd Grassmannian::Evaluate( double t ) {
	if( t == 0.0 ) return Point;

	// Step along the geodesic
	svType DH = svd.singularValues() * t;

	// Update
	MatrixXd Wt = Point;
	Wt *= svd.matrixV() * DH.array().cos().matrix().asDiagonal(); 
	Wt += svd.matrixU() * DH.array().sin().matrix().asDiagonal();
	Wt *= svd.matrixV().transpose();

	return Wt;
}

MatrixXd Grassmannian::ParallelTranslate( const MatrixXd &V, double t ) {
	if( t == 0.0 ) return V;
	
	// Step along the geodesic
	svType DH = svd.singularValues() * t;
	MatrixXd UtV = svd.matrixU().transpose() * V;

	// Update
	MatrixXd Vt = Point;
	Vt *= -svd.matrixV() * DH.array().sin().matrix().asDiagonal(); 
	Vt += svd.matrixU() * DH.array().cos().matrix().asDiagonal();
	Vt *= UtV;
	Vt += V;
	Vt -= svd.matrixU() * UtV;

	return Vt;
}

MatrixXd Grassmannian::HorizontalComponent( const MatrixXd& W, const MatrixXd& V ) {
	return V - W * (W.adjoint() * V);
}

MatrixXd Grassmannian::VerticalComponent( const MatrixXd& W, const MatrixXd& V ) {
	return W * (W.adjoint() * V);
}

