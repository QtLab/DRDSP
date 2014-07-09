#ifndef INCLUDED_DYNAMICS_RBF_FAMILY_PRODUCER
#define INCLUDED_DYNAMICS_RBF_FAMILY_PRODUCER
#include "parameter_map_producer.h"
#include "rbf_family.h"
#include <cmath>
#include <iostream>
#include <random>
#include "../data/histogram.h"
#include "../data/aabb.h"

using namespace std;

namespace DRDSP {

	template<typename F = RBF<ThinPlateSpline>>
	struct RBFFamilyProducer : ProducerBase {
		typedef PMapFamily<RBFFamily<F>,AffineXd> ReducedFamily;
		double boxScale;
		uint32_t numRBFs;

		RBFFamilyProducer() : RBFFamilyProducer(30) {}

		explicit RBFFamilyProducer( uint32_t nRBFs ) : numRBFs(nRBFs), boxScale(1.5) {}

		ReducedFamily BruteForce( const ReducedDataSystem& data, const vector<VectorXd>& parameters, uint32_t numIterations ) const {
			double Sft = 0.0, Sf = -1.0;
			uint32_t dim = data.reducedData[0].dimension;
			AABB box = data.ComputeBoundingBox();
			box.Scale(boxScale);
			RBFFamily<F> reduced( dim, numRBFs );
			ReducedFamily best;
			AffineXd A;
			vector<double> costs(numIterations);
			mt19937 mt;
			ParameterMapProducer<RBFFamily<F>> pmp;

			for(uint32_t i=0;i<numIterations;++i) {
				SetPointsRandom( reduced.centres, box, mt );
				A = pmp.SolveOrig( reduced, data, parameters );
				Sft = ComputeTotalCost( ReducedFamily( reduced, A ), data, parameters );
				costs[i] = Sft;
				if( Sft < Sf || i==0 ) {
					Sf = Sft;
					best = ReducedFamily( reduced, A );
					cout << i << " \t" << Sf << endl;
				}
			}

			HistogramGenerator histogramGenerator;
			histogramGenerator.logScale = true;
			histogramGenerator.clampMax = 0.02;
			histogramGenerator.Generate(costs).WriteCSV("output/costs.csv");

			return best;
		}

		ReducedFamily BruteForce( const ReducedDataSystem& data, const vector<VectorXd>& parameters, uint32_t numIterations, uint32_t numThreads ) const {
			vector<ReducedFamily> best(numThreads);
			AABB box = data.ComputeBoundingBox();
			box.Scale(boxScale);
			uint32_t k = data.reducedData[0].dimension * numRBFs;
			vector<future<void>> futures(numThreads);
			uint32_t iterationsPerThread = numIterations / numThreads;
			for(uint32_t i=0;i<numThreads;++i) {
				futures[i] = async( launch::async,
					[=,&data,&box,&parameters]( ReducedFamily& out ){
						out = BruteForce( data,
										  box,
										  parameters,
										  mt19937::default_seed + i * k,
										  iterationsPerThread );
					}, ref(best[i])
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

		ReducedFamily BruteForce( const ReducedDataSystem& data, const AABB& box, const vector<VectorXd>& parameters, uint32_t seed, uint32_t numIterations ) const {
			double cost = 0.0, bestCost = -1.0;
			uint32_t dim = data.reducedData[0].dimension;
			
			RBFFamily<F> reduced( dim, numRBFs );
			ReducedFamily best;
			AffineXd A;
			mt19937 mt(seed);
			ParameterMapProducer<RBFFamily<F>> pmp;

			for(uint32_t i=0;i<numIterations;++i) {
				SetPointsRandom( reduced.centres, box, mt );
				A = pmp.SolveOrig( reduced, data, parameters );
				cost = ComputeTotalCost( ReducedFamily( reduced, A ), data, parameters );
				if( cost < bestCost || i==0 ) {
					bestCost = cost;
					best = ReducedFamily( reduced, A );
					cout << i << " \t" << bestCost << endl;
				}
			}
			return best;
		}

	};

}

#endif

