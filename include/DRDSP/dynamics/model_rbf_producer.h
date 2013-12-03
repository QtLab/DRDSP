#ifndef INCLUDED_DYNAMICS_MODEL_RBF_PRODUCER
#define INCLUDED_DYNAMICS_MODEL_RBF_PRODUCER
#include "model_rbf.h"
#include "reduced_data.h"

namespace DRDSP {

	struct ModelRBFProducer {
		double fitWeight[2];
		uint16_t numRBFs;
	
		ModelRBFProducer();
		ModelRBFProducer( uint16_t nRBFs );
		double ComputeTotalCost( ModelRBF& model, const ReducedData& data ) const;

		ModelRBF ComputeModelRBF( const ReducedData& data );
		void Fit( ModelRBF& model, const ReducedData& data ) const;
		ModelRBF BruteForce( const ReducedData& data, uint32_t numIterations ) const;

	};
/*
	struct ReducedModel : ModelParameterized {
		void Randomize( const ReducedData& data );
	};

	template<typename R>
	struct ModelProducer {
		double fitWeight[2];
	
		ModelProducer();
		double ComputeTotalCost( R& model, const ReducedData& data ) const;

		R ComputeModel( const ReducedData& data ) const;
		void Fit( R& model, const ReducedData& data ) const;
		R BruteForce( const ReducedData& data, uint32_t numIterations ) const;

	};
	*/
}

#endif

