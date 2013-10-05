#include <cmath>
#include <iostream>
#include <sstream>
#include <DRDSP/dynamics/model_reduced.h>
#include <DRDSP/dynamics/radial_basis.h>
#include <Eigen/LU>
#include <DRDSP/geometry/metric.h>

using namespace std;
using namespace DRDSP;

ModelReduced::ModelReduced() : dimension(0), parameterDimension(0), numRBFs(0), rbfs(nullptr) {
}

ModelReduced::~ModelReduced() {
}

void ModelReduced::Create( uint32_t dim, uint32_t paramDim, uint32_t nRBFs ) {
	dimension = dim;
	parameterDimension = paramDim;
	numRBFs = nRBFs;
	rbfs = new RadialFunction [numRBFs];
	affine.Init(dimension,numRBFs,parameterDimension);
}

void ModelReduced::Destroy() {
	delete[] rbfs;
	rbfs = nullptr;
}


VectorXd ModelReduced::Evaluate( const VectorXd &x, const VectorXd &parameter ) const {
	MatrixXd z = affine.Evaluate(parameter);
	
	ModelRBF model( dimension, numRBFs );
	model.linear = z.block(0,0,dimension,dimension);

	for(uint32_t i=0;i<numRBFs;i++) {
		model.weights[i] = z.col(dimension+i);
		model.rbfs[i].centre = rbfs[i].centre;
	}
	return model.VectorField(x);
}

void ModelReduced::OutputText( const char *filename ) const {
	
	ofstream out;
	stringstream outfn;

	// Output results
	outfn.str("");
	outfn << "output/" << filename << "-reducedModel.csv";
	//cout << outfn.str().c_str() << endl;
	out.open(outfn.str().c_str());
	out.precision(16);
	out << dimension << " " << numRBFs << " " << parameterDimension << endl;
	for(int i=0;i<affine.coeffs.rows();i++) {	
		for(int j=0;j<affine.coeffs.cols();j++)
			out << affine.coeffs(i,j) << " ";
		out << endl;
	}
	for(uint k=0;k<numRBFs;k++) {
		for(uint j=0;j<dimension;j++)
			out << rbfs[k].centre(j) << " ";
		out << endl;
	}
	out.close();
	
}
