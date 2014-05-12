#ifndef INCLUDED_DYNAMICS_REDUCED_DATA_SYSTEM
#define INCLUDED_DYNAMICS_REDUCED_DATA_SYSTEM
#include <future>
#include "reduced_data.h"
#include "../data/data_system.h"

using namespace std;

namespace DRDSP {

	struct ReducedDataSystem {
		vector<ReducedData> reducedData;
		uint32_t numParameters;

		ReducedDataSystem();
		ReducedDataSystem( uint32_t N );
		void Create( uint32_t N );

		template<typename Family>
		ReducedDataSystem& ComputeData( Family&& family, const DataSystem& data, const MatrixXd& W ) {
			Create( data.numParameters );
			for(uint32_t i=0;i<numParameters;++i) {
				reducedData[i].ComputeData( family(data.parameters[i]), data.dataSets[i], W );
			}
			return *this;
		}

		template<typename Family>
		ReducedDataSystem& ComputeData( Family&& family, const DataSystem& data, const MatrixXd& W, uint32_t numThreads ) {
			Create( data.numParameters );
			
			vector<future<void>> futures(numThreads);
			
			for(uint32_t i=0;i<numParameters;i+=numThreads) {
				uint32_t N = min( numParameters - i, numThreads );
				for(uint32_t j=0;j<N;++j) {
					futures[j] = async( launch::async,
						[&]( ReducedData& rData, const VectorXd& parameter, const DataSet& dataSet ) {
							rData.ComputeData( family(parameter), dataSet, W );
						},
						ref(reducedData[i+j]), cref(data.parameters[i+j]), cref(data.dataSets[i+j])
					);
				}
				for(uint32_t j=0;j<N;++j) {
					futures[j].wait();
				}
			}
			return *this;
		}

		template<typename Family>
		ReducedDataSystem& ComputeDataEmbedded( Family&& family, const DataSystem& data, const MatrixXd& W ) {
			Create( data.numParameters );
			for(uint32_t i=0;i<numParameters;++i) {
				reducedData[i].ComputeDataEmbedded( family(data.parameters[i]), data.dataSets[i], W );
			}
			return *this;
		}

		template<typename Family>
		ReducedDataSystem& ComputeDataEmbedded( Family&& family, const DataSystem& data, const MatrixXd& W, uint32_t numThreads ) {
			Create( data.numParameters );
			
			vector<future<void>> futures(numThreads);
			
			for(uint32_t i=0;i<numParameters;i+=numThreads) {
				uint32_t N = min( numParameters - i, numThreads );
				for(uint32_t j=0;j<N;++j) {
					futures[j] = async( launch::async,
						[&]( ReducedData& rData, const VectorXd& parameter, const DataSet& dataSet ) {
							rData.ComputeDataEmbedded( family(parameter), dataSet, W );
						},
						ref(reducedData[i+j]), cref(data.parameters[i+j]), cref(data.dataSets[i+j])
					);
				}
				for(uint32_t j=0;j<N;++j) {
					futures[j].wait();
				}
			}
			return *this;
		}

		AABB ComputeBoundingBox() const;
		void WritePointsCSV( const char* filePrefix, const char* fileSuffix ) const;
		void WriteVectorsCSV( const char* filePrefix, const char* fileSuffix ) const;
		void WriteDerivativesCSV( const char* filePrefix, const char* fileSuffix ) const;
	};

	void Compare( const ReducedDataSystem& reducedData, const DataSystem& rdata );

}

#endif
