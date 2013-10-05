#include <cmath>
#include <iostream>
#include <sstream>
#include <Eigen/LU>
#include <DRDSP/dynamics/model_reduced_producer.h>
#include <DRDSP/dynamics/model_reduced.h>
#include <DRDSP/dynamics/affineParameterMap.h>

using namespace std;
using namespace DRDSP;

void ModelReducedProducer::GetAffineMatrices( const ReducedDataSystem& data, uint32_t numRBFs ) {

	MatrixXd y1, y2, A1, A2, Lambda, X, Y;
	uint32_t dimension = data.reducedData[0].dimension;
	uint32_t m = dimension + numRBFs;
	VectorXd temp;
	double rnorm;

	ModelRBF model( dimension, numRBFs );

	A1.setZero(m,data.count);
	y1.setZero(dimension,data.count);
	y2.setZero(dimension,data.count*dimension);
	A2.setIdentity(m,data.count*dimension);

	for(uint32_t j=0;j<data.count;j++) {
		for(uint32_t i=0;i<dimension;i++) {
			A1(i,j) = data.points[j](i);
			y1(i,j) = data.vectors[j](i);
			for(uint32_t k=0;k<dimension;k++)
				y2(i,dimension*j+k) = data.derivatives[j](i,k);
		}
		for(uint32_t i=0;i<numRBFs;i++) {
			A1(i+dimension,j) = model.rbfs[i](data.points[j]);
			temp = model.rbfs[i].Derivative(data.points[j]);
			for(uint32_t k=0;k<dimension;k++)
				A2(i+dimension,dimension*j+k) = temp(k);
		}
	}

	AffineParameterMap P;
	P.Init(dimension,m,original->pDim);
	Lambda = P.GetLambda(data->param);

	Y = A1 * y1.transpose() * (fitWeight[0]/costScale[0]) + A2 * y2.transpose() * (fitWeight[1]/costScale[1]);
	X = A1 * A1.transpose() * (fitWeight[0]/costScale[0]) + A2 * A2.transpose() * (fitWeight[1]/costScale[1]);

	B = Lambda * Y;
	A = Lambda * X * Lambda.transpose();
	
}

bool ModelReducedProducer::FitAffine() {
	MatrixXd A, B, At, Bt, z;

	uint32_t m = dimension + numRBFs;

	A.setZero(m*(original->pDim+1),m*(original->pDim+1));
	B.setZero(m*(original->pDim+1),dimension);

	for(uint32_t i=0;i<numParmameters;i++) {
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

	for(uint32_t i=0;i<numParmameters;i++) {
		z = affine.Evaluate(intermediate[i]->data->param);
		intermediate[i]->L = z.block(0,0,dim,dim);
		for(uint32_t j=0;j<numRBFs;j++) {
			intermediate[i]->weight[j] = z.col(dim+j);
		}
		//intermediate[i]->VectorizeParms();
	}
	return true;
}
