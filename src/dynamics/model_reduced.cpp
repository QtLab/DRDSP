#include <cmath>
#include <iostream>
#include <sstream>
#include <DRDSP/dynamics/model_reduced.h>
#include <DRDSP/dynamics/radial_basis.h>
#include <Eigen/LU>
#include <DRDSP/geometry/metric.h>

using namespace std;
using namespace DRDSP;

ModelReduced::ModelReduced() {
}

ModelReduced::~ModelReduced() {
}

VectorXd ModelReduced::Evaluate( const VectorXd &x, const VectorXd &b ) const {
	MatrixXd z = affine.Evaluate(b);
	
	ModelRBF model( dimension, numRBFs );
	model.linear = z.block(0,0,dimension,dimension);

	for(uint32_t i=0;i<numRBFs;i++) {
		model.weights[i] = z.col(dimension+i);
		model.rbfs[i].centre = z.col(dimension+numRBFs+i);
	}
	return model.VectorField(x);
}

void ModelReduced::GetScalesData() {
	for(uint i=0;i<numParms;i++) {
		intermediate[i]->GetScales();
		cout << "Scale factors: \t" << intermediate[i]->costScale[0] << " \t" << intermediate[i]->costScale[1] << endl;
	}
	cout << endl;
}



bool ModelReduced::FitAffine() {
	MatrixXd A, B, At, Bt, z;

	uint32_t m = dim + numRadialBasis;

	A.setZero(m*(original->pDim+1),m*(original->pDim+1));
	B.setZero(m*(original->pDim+1),dim);

	for(uint32_t i=0;i<numParms;i++) {
		intermediate[i]->GetAffineMatrices();
		B += intermediate[i]->B;
		A += intermediate[i]->A;
	}

	Eigen::FullPivLU<MatrixXd> lu(A);
	if( !lu.isInjective() ) {
		cout << "Matrix not injective, rank = " << lu.rank() << " != (" << lu.matrixLU().rows() << "," << lu.matrixLU().cols() << ")" << endl;
		return false;
	} else {
		affine.coeffs = lu.solve(B).transpose();
	}

	for(uint32_t i=0;i<numParms;i++) {
		z = affine.Evaluate(intermediate[i]->data->param);
		intermediate[i]->L = z.block(0,0,dim,dim);
		for(uint32_t j=0;j<numRadialBasis;j++) {
			intermediate[i]->weight[j] = z.col(dim+j);
		}
		//intermediate[i]->VectorizeParms();
	}
	return true;
}

void ModelReduced::LoadCentres(char *name) {
	ifstream in(name);
	if( !in ) return;

	for(uint k=0;k<numRadialBasis;k++)
		for(uint j=0;j<dim;j++)
			in >> intermediate[0]->centre[k](j);
	in.close();

	for(uint i=1;i<numParms;i++)
		for(uint k=0;k<numRadialBasis;k++)
			intermediate[i]->centre[k] = intermediate[0]->centre[k];

}

void ModelReduced::Output(char *name) {
	
	ofstream out;
	stringstream outfn;

	// Output results
	outfn.str("");
	outfn << "../../../Output/" << name << "-reducedModel.csv";
	//cout << outfn.str().c_str() << endl;
	out.open(outfn.str().c_str());
	out.precision(16);
	out << dim << " " << numRadialBasis << " " << original->pDim << endl;
	for(int i=0;i<affine.coeffs.rows();i++) {	
		for(int j=0;j<affine.coeffs.cols();j++)
			out << affine.coeffs(i,j) << " ";
		out << endl;
	}
	for(uint k=0;k<numRadialBasis;k++) {
		for(uint j=0;j<dim;j++)
			out << intermediate[0]->centre[k](j) << " ";
		out << endl;
	}
	out.close();

	outfn.str("");
	outfn << "../../../Output/" << name << "-reducedModel2.csv";
	//cout << outfn.str().c_str() << endl;
	out.open(outfn.str().c_str());
	out.precision(16);
	out << dim << " " << numRadialBasis << " " << original->pDim << " " << numParms << endl;
	for(uint i=0;i<numParms;i++) {
		for(uint j=0;j<original->pDim;j++)
			out << intermediate[i]->data->param(j) << " ";
		out << endl;
	}
	for(uint i=0;i<numParms;i++) {
		for(uint k=0;k<dim;k++)
			for(uint j=0;j<dim;j++)
				out << intermediate[i]->L(j,k) << " ";
		for(uint k=0;k<numRadialBasis;k++)
			for(uint j=0;j<dim;j++)
				out << intermediate[i]->weight[k](j) << " ";
		for(uint k=0;k<numRadialBasis;k++)
			for(uint j=0;j<dim;j++)
				out << intermediate[i]->centre[k](j) << " ";
		out << endl;
	}
	out.close();


	// Output results
	for(uint i=0;i<proj->data->numFiles;i++) {	
		outfn.str("");
		outfn << proj->data->file[i].path << "-" << name << "-parms.csv";
		out.open(outfn.str().c_str());
		out.precision(16);
		for(uint j=0;j<dim;j++) {
			for(uint k=0;k<dim;k++)
				out << intermediate[i]->L(j,k) << ",";
			out << endl;
		}
		out << endl;
		for(uint k=0;k<numRadialBasis;k++) {
			for(uint j=0;j<dim;j++)
				out << intermediate[i]->weight[k](j) << ",";
			out << endl;
		}
		out << endl;
		for(uint k=0;k<numRadialBasis;k++) {
			for(uint j=0;j<dim;j++)
				out << intermediate[i]->centre[k](j) << ",";
			out << endl;
		}
		out.close();
	}

	outfn.str("");
	outfn << "../../../Output/" << name << "-intermediates.csv";
	//cout << outfn.str().c_str() << endl;
	out.open(outfn.str().c_str());
	out.precision(16);
	for(uint i=0;i<proj->data->numFiles;i++) {	
		for(uint k=0;k<dim;k++)
			for(uint j=0;j<dim;j++)
				out << intermediate[i]->L(j,k) << ",";
		for(uint k=0;k<numRadialBasis;k++)
			for(uint j=0;j<dim;j++)
				out << intermediate[i]->weight[k](j) << ",";
		for(uint k=0;k<numRadialBasis;k++)
			for(uint j=0;j<dim;j++)
				out << intermediate[i]->centre[k](j) << ",";
		out << endl;
	}
	out.close();



	for(uint k=0;k<numParms;k++) {	
		outfn.str("");
		outfn << proj->data->file[k].path << "-" << name << "-vf.csv";
		out.open(outfn.str().c_str());
		VectorXd temp;
		for(uint32_t i=0;i<intermediate[k]->dataN;i++) {
			for(uint32_t j=0;j<dim;j++)
				out << intermediate[k]->dataQ[i](j) << ",";
			out << ",";
			temp = intermediate[k]->Evaluate(intermediate[k]->dataP[i]);
			for(uint32_t j=0;j<dim;j++)
				out << temp(j) << ",";
			out << endl;
		}
		out.close();
	}

	for(uint k=0;k<numParms;k++) {	
		outfn.str("");
		outfn << proj->data->file[k].path << "-" << name << "-vfd.csv";
		out.open(outfn.str().c_str());
		MatrixXd temp;
		VectorXd temp2, temp3;
		for(uint32_t i=0;i<intermediate[k]->dataN;i++) {
			for(uint32_t j=0;j<dim;j++)
				for(uint32_t a=0;a<dim;a++)
					out << intermediate[k]->dataR[i](j,a) << ",";
			out << ",";
			temp = intermediate[k]->EvaluateD(intermediate[k]->dataP[i]);
			for(uint32_t j=0;j<dim;j++)
				for(uint32_t a=0;a<dim;a++)
					out << temp(j,a) << ",";
			out << ",";

			temp3 = intermediate[k]->dataR[i] * intermediate[k]->dataQ[i];
			for(uint32_t j=0;j<dim;j++)
				out << temp3(j) << ",";

			temp2 = temp * intermediate[k]->Evaluate(intermediate[k]->dataP[i]);
			out << ",";
			for(uint32_t j=0;j<dim;j++)
				out << temp2(j) << ",";
			out << endl;
		}
		out.close();
	}
	
}
