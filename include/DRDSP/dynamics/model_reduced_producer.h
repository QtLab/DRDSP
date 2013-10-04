#ifndef INCLUDED_DYNAMICS_MODEL_REDUCED_PRODUCER
#define INCLUDED_DYNAMICS_MODEL_REDUCED_PRODUCER
#include "model_orig.h"
#include "model_reduced.h"
#include "reduced_data_system.h"

namespace DRDSP {

	struct ModelReducedProducer {
		double costScale[2];
		double fitWeight[2];
	
		ModelReducedProducer();
		double GetTotalCost( const ModelReduced& model, const ReducedDataSystem& data ) const;
		void ComputeScales( const ReducedDataSystem& data );

		ModelReduced ComputeModelRBF( const ReducedDataSystem& data, uint32_t numRadialBasis );
		
	};

}

#endif

