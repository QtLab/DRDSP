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
	uint32_t m = data.dimension + numRBFs;
	VectorXd temp;
	double rnorm;

	ModelRBF model( data.dimension, numRBFs );

	A1.setZero(m,data.count);
	y1.setZero(data.dimension,data.count);
	y2.setZero(data.dimension,data.count*data.dimension);
	A2.setIdentity(m,data.count*data.dimension);

	for(uint32_t j=0;j<data.count;j++) {
		for(uint32_t i=0;i<data.dimension;i++) {
			A1(i,j) = data.points[j](i);
			y1(i,j) = data.vectors[j](i);
			for(uint32_t k=0;k<data.dimension;k++)
				y2(i,data.dimension*j+k) = data.derivatives[j](i,k);
		}
		for(uint32_t i=0;i<numRBFs;i++) {
			A1(i+data.dimension,j) = model.rbfs[i](data.points[j]);
			temp = model.rbfs[i].Derivative(data.points[j]);
			for(uint32_t k=0;k<data.dimension;k++)
				A2(i+data.dimension,data.dimension*j+k) = temp(k);
		}
	}

	AffineParameterMap P;
	P.Init(data.dimension,m,original->pDim);
	Lambda = P.GetLambda(data->param);

	Y = A1 * y1.transpose() * (fitWeight[0]/costScale[0]) + A2 * y2.transpose() * (fitWeight[1]/costScale[1]);
	X = A1 * A1.transpose() * (fitWeight[0]/costScale[0]) + A2 * A2.transpose() * (fitWeight[1]/costScale[1]);

	B = Lambda * Y;
	A = Lambda * X * Lambda.transpose();
	
}
