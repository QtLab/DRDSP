#ifndef INCLUDED_DYNAMICS_FAMILY_PRODUCER
#define INCLUDED_DYNAMICS_FAMILY_PRODUCER
#include "reduced_data_system.h"
#include "../eigen_pinv.h"

using namespace std;

namespace DRDSP {
	struct FamilyProducer {
		double fitWeight[2];

		FamilyProducer() {
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
	
	

	struct GenericFamilyProducer : FamilyProducer {

		template<typename Family>
		MatrixXd ComputeParameterMap( Family&& family, const ReducedDataSystem& data, const vector<VectorXd>& parameters, const MatrixXd& warp ) {

			MatrixXd linear, linearD, D, LHS, RHS, BD, CD;
			VectorXd trans, dV, transD;

			size_t N = data.TotalPoints();

			LHS.setZero( N * data.reducedData[0].dimension * (1 + N), family.pfdim );
			RHS.setZero( LHS.rows() );

			MatrixXd rhoInv = pseudoInverse( rho );
			int64_t k = 0;
			for(uint32_t i=0;i<data.numParameters;++i) {
				D = ComputeD( family, parameter[i] );

				double scale = fitWeight[0] / data.reducedData[i].scales[0];
				double scaleD = fitWeight[1] / data.reducedData[i].scales[1];

				for(uint32_t j=0;j<data.reducedData[i].count;++j) {
					const VectorXd& x = data.reducedData[i].points[j];
					linear = ComputeLinear( family, x );
					trans = ComputeTranslation( family, x );
					auto V = rho * data.reducedData[i].vectors[j];
					BD = linear * D;
					LHS.block(k*D.cols(),0,D.cols(),D.cols()) = (BD.transpose() * BD) * scale;
					RHS.segment(k*D.cols(),D.cols()) = (BD.transpose() * (V - trans)) * scale;
					++k;
				}
				for(uint32_t j=0;j<data.reducedData[i].count;++j) {
					const VectorXd& x = data.reducedData[i].points[j];
					linearD = ComputeLinearDerivative( family, x );
					transD = ComputeTranslationDerivative( family, x );
					dV = Vectorize( rho * data.reducedData[i].derivatives[j] * rhoInv );
					CD = linearD * D;
					LHS.block(k*D.cols(),0,D.cols(),D.cols()) = (CD.transpose() * CD) * scaleD;
					RHS.segment(k*D.cols(),D.cols()) = (CD.transpose() * (dV - transD)) * scaleD;
					++k;
				}
			}
			return FullPivLU<MatrixXd>(LHS).solve(RHS);
		}

		template<typename Family>
		MatrixXd ComputeParameterMap2( Family&& family, const ReducedDataSystem& data, const vector<VectorXd>& parameters, const MatrixXd& warp ) {

			MatrixXd linear, linearD, D, LHS, transD, BD, CD;
			VectorXd trans, V, dV, RHS;

			size_t N = data.TotalPoints();

			LHS.setZero( N * data.reducedData[0].dimension * (1 + N), family.pfdim );
			RHS.setZero( LHS.rows() );

			MatrixXd rhoInv = pseudoInverse( rho );
			int64_t k = 0;
			for(uint32_t i=0;i<data.numParameters;++i) {
				D = ComputeD( family, parameter[i] );

				double scale = sqrt( fitWeight[0] / data.reducedData[i].scales[0] );
				double scaleD = sqrt( fitWeight[1] / data.reducedData[i].scales[1] );

				for(uint32_t j=0;j<data.reducedData[i].count;++j) {
					const VectorXd& x = data.reducedData[i].points[j];
					linear = ComputeLinear( family, x );
					trans = ComputeTranslation( family, x );
					V = rho * data.reducedData[i].vectors[j];
					LHS.block(k*B.rows(),0,B.rows(),D.cols()) = linear * D * scale;
					RHS.segment(k*V.size(),V.size()) = (V - trans) * scale;
					++k;
				}
				for(uint32_t j=0;j<data.reducedData[i].count;++j) {
					const VectorXd& x = data.reducedData[i].points[j];
					linearD = ComputeLinearDerivative( family, x );
					transD = ComputeTranslationDerivative( family, x );
					dV = Vectorize(rho * data.reducedData[i].derivatives[j] * rhoInv);
					LHS.block(k*C.rows(),0,C.rows(),D.cols()) = linearD * D * scaleD;
					RHS.segment(k*dV.size(),dV.size()) = (dV - transD) * scaleD;
					++k;
				}
			}
			return JacobiSVD<MatrixXd>(LHS).solve(RHS);
		}

	protected:

		template<typename Family>
		MatrixXd ComputeWarp( Family&& family, const ReducedDataSystem& data, const vector<VectorXd>& parameters ) {
			size_t cols = 0;
			for(uint32_t i=0;i<data.numParameters;++i) {
				cols += data.reducedData[i].count;
			}
			uint32_t dim = data.reducedData[0].dimension;

			MatrixXd A(dim,dim*cols);
			MatrixXd B(dim,dim*cols);
			
			typename Family::Model model;
			size_t k = 0;
			for(uint32_t i=0;i<data.numParameters;++i) {
				model = family(parameters[i]);
				const ReducedData& r = data.reducedData[i];
				for(size_t j=0;j<r.count;++j) {
					A.block(0,k*dim,dim,dim) = r.vectors[j] * r.vectors[j].transpose();
					B.block(0,k*dim,dim,dim) = r.vectors[j] * model(r.points[j]).transpose();
					++k;
				}
			}

			return FullPivLU<MatrixXd>(A).solve(B).transpose();
		}

	};

}

#endif

