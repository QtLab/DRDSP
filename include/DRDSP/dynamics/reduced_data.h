#ifndef INCLUDED_DYNAMICS_REDUCED_DATA
#define INCLUDED_DYNAMICS_REDUCED_DATA
#include "model.h"
#include "../data/data_set.h"
#include "../data/aabb.h"

#pragma warning( disable : 4510 ) // default constructor could not be generated
#pragma warning( disable : 4610 ) // can never be instantiated - user defined constructor required

#include <Eigen/SVD>

#pragma warning( default : 4610 )
#pragma warning( default : 4510 )

namespace DRDSP {

	double eps( double x );

	struct ReducedData {
		vector<VectorXd> points, vectors;
		vector<MatrixXd> derivatives;
		double scales[2];
		size_t count = 0;
		uint32_t dimension = 0;

		ReducedData() = default;
		ReducedData( uint32_t dim, size_t numPoints );
		void Create( uint32_t dim, size_t numPoints );
		AABB ComputeBoundingBox() const;
		double ComputeVectorScale();
		double ComputeDerivativeScale();
		const ReducedData& WriteData( const char* filename ) const;
		bool ReadData( const char* filename );
		const ReducedData& WritePointsCSV( const char* filename ) const;
		const ReducedData& WriteVectorsCSV( const char* filename ) const;
		const ReducedData& WriteDerivativesCSV( const char* filename ) const;

		template<typename Model>
		ReducedData& ComputeData( Model&& original, const DataSet& data, const MatrixXd& W ) {
			Create( (uint32_t)W.cols(), data.points.size() );

			for(uint32_t i=0;i<count;++i) {
				points[i].noalias() = W.adjoint() * data.points[i];
			}

			for(uint32_t i=0;i<count;++i) {
				vectors[i].noalias() = W.adjoint() * original( data.points[i] );
			}

			for(uint32_t i=0;i<count;++i) {
				derivatives[i].noalias() = W.adjoint() * original.Partials( data.points[i] ) * W;
			}
			scales[0] = ComputeVectorScale();
			scales[1] = ComputeDerivativeScale();
			return *this;
		}

		template<typename Model,typename Embedded>
		ReducedData& ComputeDataEmbedded( const ModelEmbedded<Model,Embedded>& original, const DataSet& data, const MatrixXd& W ) {
			Create( (uint32_t)W.cols(), data.points.size() );

			static const double stabilityFactor = 1.0;

			for(uint32_t i=0;i<count;++i) {
				points[i].noalias() = W.adjoint() * original.embedding(data.points[i]);
			}

			for(uint32_t i=0;i<count;++i) {
				vectors[i].noalias() = W.adjoint() * original( data.points[i] );
			}

			MatrixXd A;
			JacobiSVD<MatrixXd> svd;
			for(uint32_t i=0;i<count;++i) {
				svd.compute( W.adjoint() * original.embedding.Derivative(data.points[i]), ComputeThinU );
				uint32_t rank = 0;
				double tolerance = original.model.stateDim * eps(svd.singularValues()(0));
				for(int j=0;j<svd.nonzeroSingularValues();++j) {
					if( svd.singularValues()(j) > tolerance )
						++rank;
					else break;
				}
				A.noalias() = svd.matrixU().leftCols(rank) * svd.matrixU().leftCols(rank).transpose();
				derivatives[i].noalias() = W.adjoint() * original.Partials( data.points[i] ) * original.embedding.DerivativeAdjoint(data.points[i]) * W * A - stabilityFactor * ( MatrixXd::Identity(A.rows(),A.cols()) - A );
			}
			scales[0] = ComputeVectorScale();
			scales[1] = ComputeDerivativeScale();
			return *this;
		}

	};

}

#endif
