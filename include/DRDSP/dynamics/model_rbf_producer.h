#ifndef INCLUDED_DYNAMICS_MODEL_RBF_PRODUCER
#define INCLUDED_DYNAMICS_MODEL_RBF_PRODUCER
#include "model_rbf.h"
#include "reduced_data.h"

namespace DRDSP {

	struct ModelRBFProducer {
		double fitWeight[2];
	
		ModelRBFProducer();
		double ComputeTotalCost( const ModelRBF& model, const ReducedData& data ) const;

		ModelRBF ComputeModelRBF( const ReducedData& data, uint16_t numRadialBasis );
		void Fit( ModelRBF& model, const ReducedData& data ) const;
		ModelRBF FitBruteForce( const ReducedData& data, uint16_t numRBFs, uint32_t numInitial ) const;

	};

}

#endif

