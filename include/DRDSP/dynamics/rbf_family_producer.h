#ifndef INCLUDED_DYNAMICS_RBF_FAMILY_PRODUCER
#define INCLUDED_DYNAMICS_RBF_FAMILY_PRODUCER
#include "model.h"
#include "rbf_family.h"
#include "reduced_data_system.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <Eigen/LU>
#include "affineParameterMap.h"
#include "../data/histogram.h"

using namespace std;

namespace DRDSP {
	template<typename F = ThinPlateSpline>
	struct RBFFamilyProducer {
		double fitWeight[2], boxScale;
		uint32_t numRBFs;

		RBFFamilyProducer() : RBFFamilyProducer(30) {}

		RBFFamilyProducer( uint32_t nRBFs ) : numRBFs(nRBFs), boxScale(1.5) {
			fitWeight[0] = 0.5;
			fitWeight[1] = 0.5;
		}

		RBFFamily<F> ComputeRBFFamily( const ReducedDataSystem& data, uint32_t parameterDimension, const vector<VectorXd>& parameters ) const {

			uint32_t dimension = data.reducedData[0].dimension;

			RBFFamily<F> reduced( dimension, parameterDimension, numRBFs );
			reduced.model.SetCentresRandom( data.ComputeBoundingBox() );

			Fit(reduced,data,parameterDimension,parameters);
	
			return reduced;
		}

		double ComputeTotalCost( const RBFFamily<F>& family, const ReducedDataSystem& data, const vector<VectorXd>& parameters ) const {
			double T = 0.0;
			for(uint32_t j=0;j<data.numParameters;++j) {
				RBFModel<F> modelRBF = family( parameters[j] );
				double S1 = 0.0;
				for(uint32_t i=0;i<data.reducedData[j].count;++i) {
					S1 += ( modelRBF(data.reducedData[j].points[i]) - data.reducedData[j].vectors[i] ).squaredNorm();
				}
				S1 /= data.reducedData[j].count;
				double S2 = 0.0;
				for(uint32_t i=0;i<data.reducedData[j].count;++i) {
					S2 += ( modelRBF.Partials(data.reducedData[j].points[i]) - data.reducedData[j].derivatives[i] ).squaredNorm();
				}
				
				S2 /= data.reducedData[j].count;
				T += (fitWeight[0]/data.reducedData[j].scales[0]) * S1 + (fitWeight[1]/data.reducedData[j].scales[1]) * S2;
			}
			return T / data.numParameters;
		}

		void Fit( RBFFamily<F>& reduced, const ReducedDataSystem& data, uint32_t parameterDimension, const vector<VectorXd>& parameters ) const {
			MatrixXd A, B, Atemp, Btemp, z;

			uint32_t dimension = data.reducedData[0].dimension;
			uint32_t m = dimension + numRBFs;

			A.setZero(m*(parameterDimension+1),m*(parameterDimension+1));
			B.setZero(m*(parameterDimension+1),dimension);

			for(uint32_t i=0;i<data.numParameters;++i) {
				ComputeMatrices( reduced.model, data.reducedData[i], parameters[i], Atemp, Btemp );
				B += Btemp;
				A += Atemp;
			}

			Eigen::FullPivLU<MatrixXd> lu(A);
			if( !lu.isInjective() ) {
				//cout << "Matrix not injective, rank = " << lu.rank() << " != (" << lu.matrixLU().rows() << "," << lu.matrixLU().cols() << ")" << endl;
			} else {
				reduced.affine.coeffs = lu.solve(B).transpose();
			}

			for(uint32_t i=0;i<data.numParameters;++i) {
				z = reduced.affine(parameters[i]);
				reduced.model.linear = z.block(0,0,dimension,dimension);
				for(uint32_t j=0;j<numRBFs;++j) {
					reduced.model.weights[j] = z.col(dimension+j);
				}
			}
		}

		RBFFamily<F> BruteForce( const ReducedDataSystem& data, uint32_t parameterDimension, const vector<VectorXd>& parameters, uint32_t numIterations ) const {
			double Sft = 0.0, Sf = -1.0;

			RBFFamily<F> reduced( data.reducedData[0].dimension, parameterDimension, numRBFs );
			RBFFamily<F> best;
			AABB box = data.ComputeBoundingBox();
			box.Scale(boxScale);
			vector<double> costs(numIterations);

			for(uint32_t i=0;i<numIterations;++i) {
				reduced.model.SetCentresRandom( box );
				Fit(reduced,data,parameterDimension,parameters);
				Sft = ComputeTotalCost(reduced,data,parameters);
				costs[i] = Sft;
				if( Sft < Sf || i==0 ) {
					Sf = Sft;
					best = reduced;
					cout << i << " \t" << Sf << endl;
				}
			}

			HistogramGenerator histogramGenerator;
			histogramGenerator.logScale = true;
			histogramGenerator.clampMax = 0.02;
			histogramGenerator.Generate(costs).WriteCSV("output/costs.csv");

			return best;
		}

	protected:

		template<typename Family>
		MatrixXd ComputeRho( Family&& family, const ReducedDataSystem& data, const vector<VectorXd>& parameters ) {
			size_t cols = 0;
			for(uint32_t i=0;i<data.numParameters;++i) {
				cols += data.reducedData[i].count;
			}
			uint32_t dim = data.reducedData[0].dimension;

			MatrixXd A(dim,dim*cols);
			MatrixXd B(dim,dim*cols);
			
			typename Family::Model model;
			size_t k = 0;
			for(uint32_t i=0;i<data.numParameters;++i) {
				model = family(parameters[i]);
				const ReducedData& r = data.reducedData[i];
				for(size_t j=0;j<r.count;++j) {
					A.block(0,k*dim,dim,dim) = r.vectors[j] * r.vectors[j].transpose();
					B.block(0,k*dim,dim,dim) = r.vectors[j] * model(r.points[j]).transpose();
					++k;
				}
			}

			FullPivLU<MatrixXd> lu(A);
			return lu.solve(B).transpose();
		}

		void ComputeMatrices( const RBFModel<F>& model, const ReducedData& data, const VectorXd& parameter, MatrixXd& A, MatrixXd& B ) const {

			MatrixXd y1, y2, A1, A2, Lambda, X, Y;
			uint32_t m = data.dimension + model.numRBFs;
			uint32_t dim = data.dimension;

			A1.setZero(m,data.count);
			y1.setZero(dim,data.count);
			y2.setZero(dim,data.count*dim);
			A2.setIdentity(m,data.count*dim);

			for(uint32_t j=0;j<data.count;++j) {
				A1.block(0,j,dim,1) = data.points[j];
			}
			for(uint32_t j=0;j<data.count;++j) {
				y1.col(j) = data.vectors[j];
			}
			for(uint32_t j=0;j<data.count;++j) {
				y2.block(0,j*dim,dim,dim) = data.derivatives[j];
			}
			for(uint32_t j=0;j<data.count;++j) {
				for(uint32_t i=0;i<model.numRBFs;++i) {
					A1(dim+i,j) = model.rbfs[i](data.points[j]);
					A2.block(dim+i,j*dim,1,dim) = model.rbfs[i].Derivative(data.points[j]).transpose();
				}
			}

			AffineParameterMap Q(data.dimension,numRBFs,(uint32_t)parameter.size());
			Lambda = Q.GetLambda(parameter);

			Y = A1 * y1.transpose() * (fitWeight[0]/data.scales[0]) + A2 * y2.transpose() * (fitWeight[1]/data.scales[1]);
			X = A1 * A1.transpose() * (fitWeight[0]/data.scales[0]) + A2 * A2.transpose() * (fitWeight[1]/data.scales[1]);

			B = Lambda * Y;
			A = Lambda * X * Lambda.transpose();
	
		}

	};

}

#endif

