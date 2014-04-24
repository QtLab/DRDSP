#ifndef INCLUDED_DYNAMICS_REDUCED_DATA
#define INCLUDED_DYNAMICS_REDUCED_DATA
#include "../types.h"
#include "model.h"
#include "../data/data_set.h"
#include "../data/aabb.h"

namespace DRDSP {

	struct ReducedData {
		vector<VectorXd> points, vectors;
		vector<MatrixXd> derivatives;
		double scales[2];
		size_t count;
		uint16_t dimension;

		ReducedData();
		ReducedData( uint16_t dim, size_t numPoints );
		void Create( uint16_t dim, size_t numPoints );
		
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
		void WritePointsCSV( const char* filename ) const;
		void WriteVectorsCSV( const char* filename ) const;
	};

}

#endif
