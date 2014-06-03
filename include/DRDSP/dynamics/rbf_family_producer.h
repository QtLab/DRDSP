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

		explicit RBFFamilyProducer( uint32_t nRBFs ) : numRBFs(nRBFs), boxScale(1.5) {
			fitWeight[0] = 0.5;
			fitWeight[1] = 0.5;
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

		RBFFamily<F> BruteForce( const ReducedDataSystem& data, const vector<VectorXd>& parameters, uint32_t parameterDimension, uint32_t numIterations ) const {
			double Sft = 0.0, Sf = -1.0;

			RBFFamily<F> reduced( data.reducedData[0].dimension, parameterDimension, numRBFs );
			RBFFamily<F> best;
			AABB box = data.ComputeBoundingBox();
			box.Scale(boxScale);
			vector<double> costs(numIterations);
			mt19937 mt;

			for(uint32_t i=0;i<numIterations;++i) {
				SetCentresRandom( reduced.model, box, mt );
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

		RBFFamily<F> BruteForce( const ReducedDataSystem& data, const vector<VectorXd>& parameters, uint32_t parameterDimension, uint32_t numIterations, uint32_t numThreads ) const {
			vector<RBFFamily<F>> best(numThreads);
			AABB box = data.ComputeBoundingBox();
			box.Scale(boxScale);
			uint32_t k = data.reducedData[0].dimension * numRBFs;
			vector<future<void>> futures(numThreads);
			uint32_t iterationsPerThread = numIterations / numThreads;
			for(uint32_t i=0;i<numThreads;++i) {
				futures[i] = async( launch::async,
					[=]( RBFFamily<F>& out, const ReducedDataSystem& dataSys, const AABB& aabb, const vector<VectorXd>& params ){
						out = BruteForce( dataSys,
										  aabb,
										  params,
										  parameterDimension,
										  mt19937::default_seed + i * k,
										  iterationsPerThread );
					}, ref(best[i]), cref(data), cref(box), cref(parameters)
				);
			}

			vector<double> costs(numThreads);

			for(uint32_t i=0;i<numThreads;++i) {
				futures[i].wait();
				costs[i] = ComputeTotalCost( best[i], data, parameters );
			}

			double bestCost = costs[0];
			uint32_t bestIndex = 0;
			for(uint32_t i=1;i<numThreads;++i) {
				if( costs[i] < bestCost ) {
					bestCost = costs[i];
					bestIndex = i;
				}
			}

			return best[bestIndex];
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

			MatrixXd y1, y2, A1, A2;
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

			AffineParameterMap Q(dim,numRBFs,(uint32_t)parameter.size());
			auto Lambda = Q.GetLambda(parameter);

			auto Y = A1 * y1.transpose() * (fitWeight[0]/data.scales[0]) + A2 * y2.transpose() * (fitWeight[1]/data.scales[1]);
			auto X = A1 * A1.transpose() * (fitWeight[0]/data.scales[0]) + A2 * A2.transpose() * (fitWeight[1]/data.scales[1]);

			B = Lambda * Y;
			A = Lambda * X * Lambda.transpose();
		}

		RBFFamily<F> BruteForce( const ReducedDataSystem& data, const AABB& box, const vector<VectorXd>& parameters, uint32_t parameterDimension, uint32_t seed, uint32_t numIterations ) const {
			double Sft = 0.0, Sf = -1.0;

			RBFFamily<F> reduced( data.reducedData[0].dimension, parameterDimension, numRBFs );
			RBFFamily<F> best;
			mt19937 mt(seed);

			for(uint32_t i=0;i<numIterations;++i) {
				SetCentresRandom( reduced.model, box, mt );
				Fit(reduced,data,parameterDimension,parameters);
				Sft = ComputeTotalCost(reduced,data,parameters);
				if( Sft < Sf || i==0 ) {
					Sf = Sft;
					best = reduced;
					cout << i << " \t" << Sf << endl;
				}
			}
			return best;
		}

	protected:

		template<typename F>
		void SetCentresRandom( RBFModel<F>& model, const AABB& box, mt19937& mt ) const {
			uniform_real_distribution<double> dist;
			VectorXd diff = box.bMax - box.bMin;
			for(uint32_t i=0;i<model.numRBFs;++i) {
				for(uint32_t j=0;j<model.dimension;++j) {
					model.rbfs[i].centre(j) = box.bMin(j) + diff(j) * dist(mt);
				}
			}
		}

	};

}

#endif

