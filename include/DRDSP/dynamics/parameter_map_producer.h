#ifndef INCLUDED_DYNAMICS_PARAMETER_MAP_PRODUCER
#define INCLUDED_DYNAMICS_PARAMETER_MAP_PRODUCER
#include "reduced_data_system.h"
#include "../eigen_affine.h"
#include "../misc.h"
#include <cmath>

using namespace std;

namespace DRDSP {

	struct ProducerBase {
		double fitWeight[2];

		ProducerBase() {
			fitWeight[0] = 0.5;
			fitWeight[1] = 0.5;
		}

		template<typename Family>
		double ComputeTotalCost( const Family& family, const ReducedDataSystem& data, const vector<VectorXd>& parameters ) const {
			double T = 0.0;
			typename Family::Model model;
			for(uint32_t j=0;j<data.numParameters;++j) {
				model = family( parameters[j] );
				
				double S1 = 0.0;
				for(uint32_t i=0;i<data.reducedData[j].count;++i) {
					S1 += ( model(data.reducedData[j].points[i]) - data.reducedData[j].vectors[i] ).squaredNorm();
				}
				S1 /= data.reducedData[j].count;
				
				double S2 = 0.0;
				for(uint32_t i=0;i<data.reducedData[j].count;++i) {
					S2 += ( model.Partials(data.reducedData[j].points[i]) - data.reducedData[j].derivatives[i] ).squaredNorm();
				}
				S2 /= data.reducedData[j].count;
				
				T += (fitWeight[0]/data.reducedData[j].scales[0]) * S1 + (fitWeight[1]/data.reducedData[j].scales[1]) * S2;
			}
			return T / data.numParameters;
		}

	};

	template<typename Family>
	struct ParameterMapProducer : ProducerBase {

		AffineXd SolveSVD( const Family& family, const ReducedDataSystem& data, const vector<VectorXd>& parameters ) {
			return VecToAffine( ComputeParameterMap( family, data, parameters ), family.paramDim );
		}

		AffineXd SolveOrig( const Family& family, const ReducedDataSystem& data, const vector<VectorXd>& parameters ) {
			MatrixXd A, Atemp;
			VectorXd B, Btemp;
			uint64_t pdim = parameters[0].size();
			uint64_t m = family.paramDim * ( pdim + 1 );

			A.setZero(m,m);
			B.setZero(m);

			for(uint32_t i=0;i<data.numParameters;++i) {
				ComputeMatrices( Atemp, Btemp, family, data.reducedData[i], parameters[i] );
				A += Atemp;
				B += Btemp;
			}

			Eigen::FullPivLU<MatrixXd> lu(A);
			if( !lu.isInjective() ) {
				//cout << "Matrix not injective, rank = " << lu.rank() << " != (" << lu.matrixLU().rows() << "," << lu.matrixLU().cols() << ")" << endl;
			}

			return VecToAffine( lu.solve(B), family.paramDim );
		}

	protected:

		void ComputeMatrices( MatrixXd& A, VectorXd& B, const Family& family, const ReducedData& data, const VectorXd& parameter ) const {
			MatrixXd A1, A2;
			VectorXd y1, y2;
			uint32_t ldim = family.paramDim;
			uint32_t dim = data.dimension;
			uint32_t dim2 = dim*dim;

			A1.setZero( data.count * dim, ldim );
			y1.setZero( data.count * dim );
			y2.setZero( data.count * dim2 );
			A2.setZero( data.count * dim2, ldim );

			for(size_t i=0;i<data.count;++i) {
				A1.block(i*dim,0,dim,ldim) = family.ComputeLinear( data.points[i] );
			}
			for(size_t i=0;i<data.count;++i) {
				y1.segment(i*dim,dim) = data.vectors[i] - family.ComputeTranslation( data.points[i] );
			}
			for(size_t i=0;i<data.count;++i) {
				A2.block(i*dim2,0,dim2,ldim) = family.ComputeLinearDerivative( data.points[i] );
			}
			for(size_t i=0;i<data.count;++i) {
				y2.segment(i*dim2,dim2) = Vectorize( data.derivatives[i] - family.ComputeTranslationDerivative( data.points[i] ) );
			}
			
			auto Y = A1.transpose() * y1 * (fitWeight[0]/data.scales[0]) + A2.transpose() * y2 * (fitWeight[1]/data.scales[1]);
			auto X = A1.transpose() * A1 * (fitWeight[0]/data.scales[0]) + A2.transpose() * A2 * (fitWeight[1]/data.scales[1]);
			auto D = ComputeD(parameter,ldim);

			B = D.transpose() * Y;
			A = D.transpose() * X * D;
		}

		VectorXd ComputeParameterMap( const Family& family, const ReducedDataSystem& data, const vector<VectorXd>& parameters ) {
			MatrixXd linear, linearD, D, LHS, BD, CD, transD;
			VectorXd trans, RHS;
			size_t N = data.TotalPoints();

			LHS.setZero(
				N * data.reducedData[0].dimension * (1 + data.reducedData[0].dimension),
				family.paramDim * ( parameters[0].size() + 1 )
			);
			RHS.setZero( LHS.rows() );

			int64_t k1 = 0, k2 = 0;
			for(uint32_t i=0;i<data.numParameters;++i) {
				D = ComputeD( parameters[i], family.paramDim );

				double scale = std::sqrt( fitWeight[0] / data.reducedData[i].scales[0] );
				double scaleD = std::sqrt( fitWeight[1] / data.reducedData[i].scales[1] );

				for(uint32_t j=0;j<data.reducedData[i].count;++j) {
					const VectorXd& x = data.reducedData[i].points[j];
					linear = family.ComputeLinear( x );
					trans = family.ComputeTranslation( x );
					const VectorXd& V = data.reducedData[i].vectors[j];
					LHS.block(k1,0,linear.rows(),D.cols()) = linear * D * scale;
					RHS.segment(k2,V.size()) = (V - trans) * scale;
					k1 += linear.rows();
					k2 += V.size();
				}
				for(uint32_t j=0;j<data.reducedData[i].count;++j) {
					const VectorXd& x = data.reducedData[i].points[j];
					linearD = family.ComputeLinearDerivative( x );
					transD = family.ComputeTranslationDerivative( x );
					const MatrixXd& dV = data.reducedData[i].derivatives[j];
					LHS.block(k1,0,linearD.rows(),D.cols()) = linearD * D * scaleD;
					RHS.segment(k2,dV.size()) = Vectorize(dV - transD) * scaleD;
					k1 += linearD.rows();
					k2 += dV.size();
				}
			}
			return JacobiSVD<MatrixXd>(LHS,ComputeThinU|ComputeThinV).solve(RHS);
		}

		static MatrixXd ComputeD( const VectorXd& parameter, int64_t ldim ) {
			MatrixXd D;
			int64_t pdim = parameter.size();
			D.setZero( ldim, ldim * ( pdim + 1 ) );
			for(int64_t i=0;i<pdim;++i) {
				D.block(0,ldim*i,ldim,ldim).setIdentity() *= parameter[i];
			}
			D.topRightCorner(ldim,ldim).setIdentity();
			return D;
		}
	};
}

#endif

