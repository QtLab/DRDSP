#include <cmath>
#include <iostream>
#include <sstream>
#include <DRDSP/dynamics/model_reduced.h>
#include <DRDSP/dynamics/radial_basis.h>
#include <Eigen/LU>
#include <DRDSP/geometry/metric.h>

using namespace std;
using namespace DRDSP;

ModelReduced::ModelReduced() : dimension(0), parameterDimension(0) {
}

ModelReduced::~ModelReduced() {
}

void ModelReduced::Create( uint16_t dim, uint8_t paramDim, uint16_t nRBFs ) {
	dimension = dim;
	parameterDimension = paramDim;
	model.Create(dim,nRBFs);
	affine.Init(dimension,nRBFs,parameterDimension);
}

void ModelReduced::Destroy() {
}

ModelRBF ModelReduced::ComputeModelRBF( const VectorXd& parameter ) {
	MatrixXd z = affine.Evaluate(parameter);

	model.linear = z.block(0,0,dimension,dimension);

	for(uint16_t i=0;i<model.numRBFs;i++) {
		model.weights[i] = z.col(dimension+i);
	}
	return model;
}

VectorXd ModelReduced::Evaluate( const VectorXd &x, const VectorXd &parameter ) {
	MatrixXd z = affine.Evaluate(parameter);

	model.linear = z.block(0,0,dimension,dimension);

	for(uint16_t i=0;i<model.numRBFs;i++) {
		model.weights[i] = z.col(dimension+i);
	}
	return model.VectorField(x);
}

void ModelReduced::OutputText( const char *filename ) const {
	
	ofstream out;
	out.open(filename);
	out.precision(16);
	out << dimension << " " << model.numRBFs << " " << parameterDimension << endl;
	for(int i=0;i<affine.coeffs.rows();i++) {	
		for(int j=0;j<affine.coeffs.cols();j++)
			out << affine.coeffs(i,j) << " ";
		out << endl;
	}
	for(uint16_t k=0;k<model.numRBFs;k++) {
		for(uint16_t j=0;j<dimension;j++)
			out << model.rbfs[k].centre(j) << " ";
		out << endl;
	}
	out.close();
	
}
