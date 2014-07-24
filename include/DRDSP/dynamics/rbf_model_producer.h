#ifndef INCLUDED_DYNAMICS_RBF_MODEL_PRODUCER
#define INCLUDED_DYNAMICS_RBF_MODEL_PRODUCER
#include "rbf_model.h"
#include "reduced_data.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <Eigen/LU>
#include "../data/histogram.h"
#include "producer_base.h"
#include "../misc.h"

namespace DRDSP {

	template<typename Family>
	struct RBFModelProducer : ProducerBase {
		typedef typename Family::Model Model;
		double boxScale = 1.5;
		uint32_t numRBFs;

		explicit RBFModelProducer( uint32_t numRBFs ) : numRBFs(numRBFs) {}

		VectorXd Solve( const Family& family, const ReducedData& data ) const {
			MatrixXd A1, A2;
			VectorXd y1, y2;
			uint32_t ldim = family.paramDim;
			uint32_t dim = data.dimension;
			uint32_t dim2 = dim*dim;

			A1.setZero( data.count * dim, ldim );
			y1.setZero( data.count * dim );
			y2.setZero( data.count * dim2 );
			A2.setZero( data.count * dim2, ldim );

			for(size_t i=0;i<data.count;++i) {
				A1.block(i*dim,0,dim,ldim) = family.ComputeLinear( data.points[i] );
			}
			for(size_t i=0;i<data.count;++i) {
				y1.segment(i*dim,dim) = data.vectors[i] - family.ComputeTranslation( data.points[i] );
			}
			for(size_t i=0;i<data.count;++i) {
				A2.block(i*dim2,0,dim2,ldim) = family.ComputeLinearDerivative( data.points[i] );
			}
			for(size_t i=0;i<data.count;++i) {
				y2.segment(i*dim2,dim2) = Vectorize( data.derivatives[i] - family.ComputeTranslationDerivative( data.points[i] ) );
			}
			
			auto A = A1.transpose() * A1 * (fitWeight[0]/data.scales[0]) + A2.transpose() * A2 * (fitWeight[1]/data.scales[1]);
			auto b = A1.transpose() * y1 * (fitWeight[0]/data.scales[0]) + A2.transpose() * y2 * (fitWeight[1]/data.scales[1]);

			Eigen::FullPivLU<MatrixXd> lu(A);
			if( !lu.isInjective() ) {
				//cout << "Matrix not injective, rank = " << lu.rank() << " != (" << lu.matrixLU().rows() << "," << lu.matrixLU().cols() << ")" << endl;
			}
			return lu.solve(b);
		}

		Model BruteForce( const ReducedData& data, uint32_t numIterations ) const {
			double cost = 0.0, bestCost = -1.0;

			Model family( data.dimension, numRBFs );
			Model model, best;
			AABB box = data.ComputeBoundingBox();
			box.Scale(boxScale);
			vector<double> costs(numIterations);
			mt19937 mt;

			for(uint32_t i=0;i<numIterations;++i) {
				SetPointsRandom( family.centres, box, mt );
				model = family( Solve( family, data ) );
				costs[i] = ComputeTotalCost( model, data );
				if( costs[i] < bestCost || i==0 ) {
					bestCost = costs[i];
					best = model;
					cout << i << " \t" << bestCost << endl;
				}
			}
			HistogramGenerator(100).Generate(costs).WriteCSV("output/costs.csv");
			return best;
		}

		Model BruteForce( const ReducedData& data, uint32_t numIterations, uint32_t numThreads ) const {
			vector<Model> best(numThreads);
			AABB box = data.ComputeBoundingBox();
			box.Scale(boxScale);
			uint32_t k = data.dimension * numRBFs;
			vector<future<void>> futures(numThreads);
			uint32_t iterationsPerThread = numIterations / numThreads;
			for(uint32_t i=0;i<numThreads;++i) {
				futures[i] = async( launch::async,
					[=,&data,&box]( Model& out ){
						out = BruteForce( data,
										  box,
										  mt19937::default_seed + i * k,
										  iterationsPerThread );
					}, ref(best[i])
				);
			}

			vector<double> costs(numThreads);

			for(uint32_t i=0;i<numThreads;++i) {
				futures[i].wait();
				costs[i] = ComputeTotalCost( best[i], data );
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

		Model BruteForce( const ReducedData& data, const AABB& box, uint32_t seed, uint32_t numIterations ) const {
			double cost = 0.0, bestCost = -1.0;
			uint32_t dim = data.dimension;
			
			Family family( dim, numRBFs );
			Model model, best;
			mt19937 mt(seed);

			for(uint32_t i=0;i<numIterations;++i) {
				SetPointsRandom( family.centres, box, mt );
				model = family( Solve( family, data ) );
				cost = ComputeTotalCost( model, data );
				if( cost < bestCost || i==0 ) {
					bestCost = cost;
					best = model;
					cout << i << " \t" << bestCost << endl;
				}
			}
			return best;
		}
	};
}

#endif

