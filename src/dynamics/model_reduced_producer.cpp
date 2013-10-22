#include <cmath>
#include <iostream>
#include <sstream>
#include <Eigen/LU>
#include <DRDSP/dynamics/model_reduced_producer.h>
#include <DRDSP/dynamics/model_reduced.h>
#include <DRDSP/dynamics/affineParameterMap.h>

using namespace std;
using namespace DRDSP;

void ModelReducedProducer::ComputeMatrices( const ModelRBF& model, const ReducedData& data, const VectorXd& parameter, MatrixXd& A, MatrixXd& B ) const {

	MatrixXd y1, y2, A1, A2, Lambda, X, Y;
	uint16_t m = data.dimension + model.numRBFs;
	VectorXd temp;

	A1.setZero(m,data.count);
	y1.setZero(data.dimension,data.count);
	y2.setZero(data.dimension,data.count*data.dimension);
	A2.setIdentity(m,data.count*data.dimension);

	for(uint32_t j=0;j<data.count;j++) {
		for(uint16_t i=0;i<data.dimension;i++) {
			A1(i,j) = data.points[j](i);
			y1(i,j) = data.vectors[j](i);
			for(uint16_t k=0;k<data.dimension;k++)
				y2(i,data.dimension*j+k) = data.derivatives[j](i,k);
		}
		for(uint16_t i=0;i<model.numRBFs;i++) {
			A1(i+data.dimension,j) = model.rbfs[i](data.points[j]);
			temp = model.rbfs[i].Derivative(data.points[j]);
			for(uint16_t k=0;k<data.dimension;k++)
				A2(i+data.dimension,data.dimension*j+k) = temp(k);
		}
	}

	AffineParameterMap P;
	P.Init(data.dimension,numRBFs,(uint8_t)parameter.size());
	Lambda = P.GetLambda(parameter);

	Y = A1 * y1.transpose() * (fitWeight[0]/data.scales[0]) + A2 * y2.transpose() * (fitWeight[1]/data.scales[1]);
	X = A1 * A1.transpose() * (fitWeight[0]/data.scales[0]) + A2 * A2.transpose() * (fitWeight[1]/data.scales[1]);

	B = Lambda * Y;
	A = Lambda * X * Lambda.transpose();
	
}

ModelReduced ModelReducedProducer::ComputeModelReduced( const ReducedDataSystem& data, uint8_t parameterDimension, const VectorXd* parameters ) const {
	MatrixXd A, B, Atemp, Btemp, z;

	uint16_t dimension = data.reducedData[0].dimension;
	uint16_t m = dimension + numRBFs;

	ModelReduced reduced;
	reduced.Create( dimension, parameterDimension, numRBFs );

	reduced.model.SetCentresRandom( data.ComputeBoundingBox() );

	A.setZero(m*(data.numParameters+1),m*(data.numParameters+1));
	B.setZero(m*(data.numParameters+1),dimension);

	for(uint16_t i=0;i<data.numParameters;i++) {
		ComputeMatrices( reduced.model, data.reducedData[i], parameters[i], Atemp, Btemp );
		B += Btemp;
		A += Atemp;
	}

	Eigen::FullPivLU<MatrixXd> lu(A);
	if( !lu.isInjective() ) {
		cout << "Matrix not injective, rank = " << lu.rank() << " != (" << lu.matrixLU().rows() << "," << lu.matrixLU().cols() << ")" << endl;
	} else {
		reduced.affine.coeffs = lu.solve(B).transpose();
	}

	for(uint16_t i=0;i<data.numParameters;i++) {
		z = reduced.affine.Evaluate(parameters[i]);
		reduced.model.linear = z.block(0,0,dimension,dimension);
		for(uint16_t j=0;j<numRBFs;j++) {
			reduced.model.weights[j] = z.col(dimension+j);
		}
	}
	return std::move(reduced);
}

double ModelReducedProducer::ComputeTotalCost( ModelReduced& model, const ReducedDataSystem& data, const VectorXd* parameters ) const {
	double T = 0.0;
	for(uint16_t j=0;j<data.numParameters;j++) {
		double S1 = 0.0, S2 = 0.0;
		ModelRBF modelRBF = model.ComputeModelRBF( parameters[j] );
		for(uint32_t i=0;i<data.reducedData[j].count;i++) {
			S1 += ( modelRBF.VectorField(data.reducedData[j].points[i]) - data.reducedData[j].vectors[i] ).squaredNorm();
			S2 += ( modelRBF.VectorFieldDerivative(data.reducedData[j].points[i]) - data.reducedData[j].derivatives[i] ).squaredNorm();
		}
		T += (fitWeight[0]/data.reducedData[j].scales[0]) * S1 + (fitWeight[1]/data.reducedData[j].scales[1]) * S2;
	}
	return T;
}

