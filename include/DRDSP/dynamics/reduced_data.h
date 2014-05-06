#ifndef INCLUDED_DYNAMICS_REDUCED_DATA
#define INCLUDED_DYNAMICS_REDUCED_DATA
#include <Eigen/SVD>
#include "model.h"
#include "../data/data_set.h"
#include "../data/aabb.h"

namespace DRDSP {

	double eps( double x );

	struct ReducedData {
		vector<VectorXd> points, vectors;
		vector<MatrixXd> derivatives;
		double scales[2];
		size_t count;
		uint32_t dimension;

		ReducedData();
		ReducedData( uint32_t dim, size_t numPoints );
		void Create( uint32_t dim, size_t numPoints );
		AABB ComputeBoundingBox() const;
		double ComputeVectorScale();
		double ComputeDerivativeScale();
		void WriteData( const char* filename ) const;
		bool ReadData( const char* filename );
		void WritePointsCSV( const char* filename ) const;
		void WriteVectorsCSV( const char* filename ) const;

		template<typename Model>
		void ComputeData( Model&& original, const DataSet& data, const MatrixXd& W ) {
			Create( (uint32_t)W.cols(), data.points.size() );

			for(uint32_t i=0;i<count;++i) {
				points[i] = W.adjoint() * data.points[i];
			}

			for(uint32_t i=0;i<count;++i) {
				vectors[i] = W.adjoint() * original( data.points[i] );
			}

			for(uint32_t i=0;i<count;++i) {
				derivatives[i] = W.adjoint() * original.Partials( data.points[i] ) * W;
			}
			scales[0] = ComputeVectorScale();
			scales[1] = ComputeDerivativeScale();
		}

		template<typename Model,typename Embedded>
		void ComputeDataEmbedded( const ModelEmbedded<Model,Embedded>& original, const DataSet& data, const MatrixXd& W ) {
			Create( (uint32_t)W.cols(), data.points.size() );

			static double stabilityFactor = 1.0;

			MatrixXd A;

			for(uint32_t i=0;i<count;++i) {
				points[i] = W.adjoint() * original.embedding(data.points[i]);
			}
			for(uint32_t i=0;i<count;++i) {
				vectors[i] = W.adjoint() * original( data.points[i] );
			}
			for(uint32_t i=0;i<count;++i) {
				JacobiSVD<MatrixXd> svd(W.adjoint() * original.embedding.Derivative(data.points[i]),ComputeThinU);
				uint32_t rank = 0;
				double tolerance = original.model.dimension * eps(svd.singularValues()(0));
				for(int j=0;j<svd.nonzeroSingularValues();++j) {
					if( svd.singularValues()(j) > tolerance )
						++rank;
					else break;
				}
				A = svd.matrixU().leftCols(rank) * svd.matrixU().leftCols(rank).transpose();
				derivatives[i] = W.adjoint() * original.Partials( data.points[i] ) * original.embedding.DerivativeAdjoint(data.points[i]) * W * A - stabilityFactor * ( MatrixXd::Identity(A.rows(),A.cols()) - A );
			}
			scales[0] = ComputeVectorScale();
			scales[1] = ComputeDerivativeScale();
		}

	};

}

#endif
