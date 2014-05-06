#ifndef INCLUDED_DYNAMICS_RBF_MODEL_PRODUCER
#define INCLUDED_DYNAMICS_RBF_MODEL_PRODUCER
#include "rbf_model.h"
#include "reduced_data.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <Eigen/LU>
#include "../data/histogram.h"

namespace DRDSP {

	template<typename F>
	struct RBFModelProducer {
		double fitWeight[2];
		uint32_t numRBFs;

		RBFModelProducer() : RBFModelProducer(30) {}

		RBFModelProducer( uint32_t nRBFs ) : numRBFs(nRBFs) {
			fitWeight[0] = 0.5;
			fitWeight[1] = 0.5;
		}

		double ComputeTotalCost( const RBFModel<F>& model, const ReducedData& data ) const {
			double S1 = 0.0;
			for(uint32_t i=0;i<data.count;++i) {
				S1 += ( model(data.points[i]) - data.vectors[i] ).squaredNorm();
			}
			S1 /= data.count;
			double S2 = 0.0;
			for(uint32_t i=0;i<data.count;++i) {
				S2 += ( model.Partials(data.points[i]) - data.derivatives[i] ).squaredNorm();
			}
			S2 /= data.count;
			return (fitWeight[0]/data.scales[0]) * S1 + (fitWeight[1]/data.scales[1]) * S2;
		}

		RBFModel<F> ComputeModelRBF( const ReducedData& data ) {

			RBFModel<F> model(data.dimension,numRBFs);
			model.SetCentresRandom( data.ComputeBoundingBox() );

			Fit(model,data);

			return model;
		}

		void Fit( RBFModel<F>& model, const ReducedData& data ) const {
			MatrixXd y1, y2, A1, A2, z1, z2, z, P, PI;
			uint32_t m = data.dimension + model.numRBFs;
			VectorXd temp;

			A1.setZero(m,data.count);
			y1.setZero(data.dimension,data.count);
			y2.setZero(data.dimension,data.count*data.dimension);
			A2.setIdentity(m,data.count*data.dimension);

			for(uint32_t j=0;j<data.count;++j) {
				for(uint32_t i=0;i<data.dimension;++i) {
					A1(i,j) = data.points[j](i);
					y1(i,j) = data.vectors[j](i);
					for(uint32_t k=0;k<data.dimension;++k)
						y2(i,data.dimension*j+k) = data.derivatives[j](i,k);
				}
				for(uint32_t i=0;i<model.numRBFs;++i) {
					A1(i+data.dimension,j) = model.rbfs[i](data.points[j]);
					temp = model.rbfs[i].Derivative(data.points[j]);
					for(uint32_t k=0;k<data.dimension;++k)
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
			for(uint32_t a=0;a<model.numRBFs;++a)
				model.weights[a] = z.col(data.dimension+a);
		}

		RBFModel<F> BruteForce( const ReducedData& data, uint32_t numIterations ) const {
			double Sft = 0.0, Sf = -1.0;

			RBFModel<F> model( data.dimension, numRBFs ), best;
			AABB box = data.ComputeBoundingBox();
			vector<double> costs(numIterations);

			for(uint32_t i=0;i<numIterations;++i) {
				model.SetCentresRandom( box );
				Fit(model,data);
				Sft = ComputeTotalCost(model,data);
				costs[i] = Sft;
				if( Sft < Sf || i==0 ) {
					Sf = Sft;
					best = model;
					cout << i << " \t" << Sf << endl;
				}
			}

			HistogramGenerator(100).Generate(costs).WriteCSV("output/costs.csv");

			return best;
		}

	};

}

#endif

