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

			for(uint32_t i=0;i<count;i++) {
				points[i] = W.adjoint() * data.points[i];
				vectors[i] = W.adjoint() * original( data.points[i] );
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

			for(uint32_t i=0;i<count;i++) {

				points[i] = W.adjoint() * original.embedding(data.points[i]);
				vectors[i] = W.adjoint() * original( data.points[i] );

				JacobiSVD<MatrixXd> svd(W.adjoint() * original.embedding.Derivative(data.points[i]),ComputeThinU);
				uint32_t rank = 0;
				double tolerance = original.model.dimension * eps(svd.singularValues()(0));
				for(int j=0;j<svd.nonzeroSingularValues();j++) {
					if( svd.singularValues()(j) > tolerance )
						rank++;
					else break;
				}
				A = svd.matrixU().leftCols(rank) * svd.matrixU().leftCols(rank).transpose();
				derivatives[i] = W.adjoint() * original.Partials( data.points[i] ) * original.embedding.DerivativeAdjoint(data.points[i]) * W * A - stabilityFactor * ( MatrixXd::Identity(A.rows(),A.cols()) - A );
			}
			scales[0] = ComputeVectorScale();
			scales[1] = ComputeDerivativeScale();
		}

		template<typename Model>
		void ComputeDataCW( Model&& original, const DataSet& data, const MatrixXd& W ) {
			Create( (uint32_t)W.cols(), data.points.size() );

			for(uint32_t i=0;i<count;i++) {

				points[i] = W.adjoint() * data.points[i];

				for(uint32_t j=0;j<dimension;j++)
					for(uint32_t k=0;k<original.dimension;k++)
						vectors[i](j) = W(k,j) * original( data.points[i], k );

				for(uint32_t c1=0;c1<original.dimension;c1++)
					for(uint32_t c2=0;c2<original.dimension;c2++) {
						double temp = original.Partials( data.points[i], c1, c2 );
						if( temp != 0.0 ) {
							for(uint32_t k=0;k<dimension;k++)
								for(uint32_t j=0;j<dimension;j++)
									derivatives[i](j,k) += W(c1,j) * temp * W(c2,k);
						}
					}
			}
			scales[0] = ComputeVectorScale();
			scales[1] = ComputeDerivativeScale();
		}

		template<typename Model,typename Embedded>
		void ComputeDataEmbeddedCW( const ModelEmbedded<Model,Embedded>& original, const DataSet& data, const MatrixXd& W ) {
			Create( (uint16_t)W.cols(), data.points.size() );

			static double stabilityFactor = 1.0;

			MatrixXd A, dPhi, tempMatrix;
	
			for(uint32_t i=0;i<count;i++) {
				dPhi.setZero(dimension,original.embedding.eDim);
				tempMatrix.setZero(dimension,dimension);

				points[i] = W.adjoint() * original.embedding(data.points[i]);
				for(uint32_t j=0;j<dimension;j++)
					for(uint32_t k=0;k<original.embedding.eDim;k++) {
						vectors[i](j) = W(k,j) * original( data.points[i], k );
						for(uint32_t c=0;c<original.embedding.oDim;c++)
							dPhi(j,c) += W(k,j) * original.embedding.Derivative(data.points[i],k,c);
					}

				JacobiSVD<MatrixXd> svd(dPhi,ComputeThinU);
				uint32_t rank = 0;
				double tolerance = original.model.dimension * eps(svd.singularValues()(0));
				for(int j=0;j<svd.nonzeroSingularValues();j++) {
					if( svd.singularValues()(j) > tolerance )
						rank++;
					else break;
				}
				A = svd.matrixU().leftCols(rank) * svd.matrixU().leftCols(rank).transpose();
		
				for(uint32_t c1=0;c1<original.model.dimension;c1++)
					for(uint32_t c2=0;c2<original.model.dimension;c2++) {
						double temp1 = original.Partials( data.points[i], c1, c2 );
						if( temp1 == 0.0 ) continue;
						for(uint32_t c3=0;c3<original.embedding.eDim;c3++) {
							double temp2 = original.embedding.DerivativeAdjoint(data.points[i],c2,c3);
							if( temp2 == 0.0 ) continue;
							double temp = temp1 * temp2;
							for(uint32_t k=0;k<dimension;k++)
								for(uint32_t j=0;j<dimension;j++)
									tempMatrix(j,k) += W(c1,j) * temp * W(c2,k);
						}
					}
				derivatives[i] = tempMatrix * A - stabilityFactor * ( MatrixXd::Identity(A.rows(),A.cols()) - A );
			}
			scales[0] = ComputeVectorScale();
			scales[1] = ComputeDerivativeScale();
		}

	};

}

#endif
