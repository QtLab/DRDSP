#ifndef INCLUDED_DYNAMICS_PRODUCER_BASE
#define INCLUDED_DYNAMICS_PRODUCER_BASE
#include "reduced_data_system.h"

namespace DRDSP {

	struct ProducerBase {
		double fitWeight[2];

		ProducerBase() {
			fitWeight[0] = 0.5;
			fitWeight[1] = 0.5;
		}

		template<typename Model>
		double ComputeTotalCost( const Model& model, const ReducedData& data ) const {
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

		template<typename Family>
		double ComputeTotalCost( const Family& family, const ReducedDataSystem& data, const vector<VectorXd>& parameters ) const {
			double T = 0.0;
			for(uint32_t j=0;j<data.numParameters;++j) {
				T += ComputeTotalCost( family( parameters[j] ), data[j] );
			}
			return T / data.numParameters;
		}
	};
}

#endif

