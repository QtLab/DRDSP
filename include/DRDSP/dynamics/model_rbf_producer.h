#ifndef INCLUDED_DYNAMICS_MODEL_RBF_PRODUCER
#define INCLUDED_DYNAMICS_MODEL_RBF_PRODUCER
#include "model_rbf.h"
#include "reduced_data.h"

namespace DRDSP {

	struct ModelRBFProducer {
		double fitWeight[2];
		uint16_t numRBFs;
	
		ModelRBFProducer();
		double ComputeTotalCost( const ModelRBF& model, const ReducedData& data ) const;

		ModelRBF ComputeModelRBF( const ReducedData& data );
		void Fit( ModelRBF& model, const ReducedData& data ) const;
		ModelRBF BruteForce( const ReducedData& data, uint32_t numIterations ) const;

	};

}

#endif

