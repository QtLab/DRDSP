#ifndef INCLUDED_DYNAMICS_REDUCED_DATA
#define INCLUDED_DYNAMICS_REDUCED_DATA
#include "../types.h"
#include "model_orig.h"
#include "../data/data_set.h"
#include "../data/aabb.h"

namespace DRDSP {

	struct ReducedData {
		VectorXd *points, *vectors;
		MatrixXd *derivatives;
		double scales[2];
		uint32_t count;
		uint16_t dimension;

		ReducedData();
		ReducedData( uint16_t dim, uint32_t numPoints );
		ReducedData( const ReducedData& rhs );
		ReducedData( ReducedData&& rhs );
		~ReducedData();
		ReducedData& operator=( const ReducedData& rhs );
		ReducedData& operator=( ReducedData&& rhs );
		void Create( uint16_t dim, uint32_t numPoints );
		void Destroy();
		
		void ComputeData( Model& model, const DataSet& data, const MatrixXd& W );
		void ComputeData( ModelCW& model, const DataSet& data, const MatrixXd& W );
		void ComputeData( ModelParameterized& model, const VectorXd& parameter, const DataSet& data, const MatrixXd& W );
		void ComputeData( ModelParameterizedCW& model, const VectorXd& parameter, const DataSet& data, const MatrixXd& W );
		void ComputeData( ModelEmbedded& model, const DataSet& data, const MatrixXd& W );
		void ComputeData( ModelParameterizedEmbedded& model, const VectorXd& parameter, const DataSet& data, const MatrixXd& W );
		void ComputeData( ModelParameterizedEmbeddedCW& model, const VectorXd& parameter, const DataSet& data, const MatrixXd& W );
		
		AABB ComputeBoundingBox() const;
		double ComputeVectorScale();
		double ComputeDerivativeScale();
		void WriteData( const char* filename ) const;
		bool ReadData( const char* filename );
		void WritePointsText( const char* filename ) const;
		void WriteVectorsText( const char* filename ) const;
	};

}

#endif
