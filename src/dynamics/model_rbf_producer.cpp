#include <cmath>
#include <iostream>
#include <sstream>
#include <Eigen/LU>
#include <DRDSP/dynamics/model_rbf_producer.h>
#include <DRDSP/dynamics/model_rbf.h>

using namespace std;
using namespace DRDSP;

ModelRBFProducer::ModelRBFProducer() {
	fitWeight[0] = 1.0;
	fitWeight[1] = 1.0;
}

double ModelRBFProducer::ComputeTotalCost( const ModelRBF& model, const ReducedData& data ) const {
	double S1 = 0.0, S2 = 0.0;
	for(uint32_t i=0;i<data.count;i++) {
		S1 += ( model.VectorField(data.points[i]) - data.vectors[i] ).squaredNorm();
		S2 += ( model.VectorFieldDerivative(data.points[i]) - data.derivatives[i] ).squaredNorm();
	}
	return (fitWeight[0]/data.scales[0]) * S1 + (fitWeight[1]/data.scales[1]) * S2;
}

ModelRBF ModelRBFProducer::ComputeModelRBF( const ReducedData& data, uint16_t numRBFs ) {

	ModelRBF model(data.dimension,numRBFs);
	model.SetCentresRandom( data.ComputeBoundingBox() );

	Fit(model,data);

	return std::move(model);
}

void ModelRBFProducer::Fit( ModelRBF& model, const ReducedData& data ) const {
	MatrixXd y1, y2, A1, A2, z1, z2, z, P, PI;
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

	P = (A1 * A1.transpose()) * (fitWeight[0]/data.scales[0]) + (A2 * A2.transpose()) * (fitWeight[1]/data.scales[1]);
	
	FullPivLU<MatrixXd> lu(P);
	PI = lu.inverse();

	z1 = y1 * A1.transpose() * PI;
	z2 = y2 * A2.transpose() * PI;
	z = z1 * (fitWeight[0]/data.scales[0]) + z2 * (fitWeight[1]/data.scales[1]);

	model.linear = z.block(0,0,data.dimension,data.dimension);
	for(uint16_t a=0;a<model.numRBFs;a++)
		model.weights[a] = z.col(data.dimension+a);
}

ModelRBF ModelRBFProducer::FitBruteForce( const ReducedData& data, uint16_t numRBFs, uint32_t numIterations ) const {
	double Sft = 0.0, Sf = -1.0;

	ModelRBF model( data.dimension, numRBFs );
	AABB box = data.ComputeBoundingBox();
	ModelRBF best;

	for(uint32_t i=0;i<numIterations;i++) {
		model.SetCentresRandom( box );
		Fit(model,data);
		Sft = ComputeTotalCost(model,data);

		if( Sft < Sf || i==0 ) {
			Sf = Sft;
			best = model;
			cout << Sf << endl;
		}
	}
	return std::move(best);
}
