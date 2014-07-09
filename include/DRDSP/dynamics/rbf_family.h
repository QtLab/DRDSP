#ifndef INCLUDED_DYNAMICS_RBF_FAMILY
#define INCLUDED_DYNAMICS_RBF_FAMILY
#include <iostream>
#include <fstream>
#include "rbf_model.h"

using namespace std;

namespace DRDSP {

	template<typename F = ThinPlateSpline>
	struct RBFNLFamily : Family<RBFModel<F>> {
		uint32_t nRBFs;

		RBFNLFamily() : nRBFs(0) {}

		RBFNLFamily( uint32_t dim, uint32_t nRBFs ) : Family<RBFModel<F>>(dim,dim*(dim+2*nRBFs)), nRBFs(nRBFs) {}

		RBFModel<F> operator()( const VectorXd& parameter ) const {
			RBFModel<F> model(stateDim,nRBFs);
			for(uint32_t i=0;i<stateDim;++i)
				model.linear.col(i) = parameter.segment(i*stateDim,stateDim);
			for(uint32_t i=0;i<nRBFs;++i)
				model.rbfs[i].weight = parameter.segment((i+stateDim)*stateDim,stateDim);
			for(uint32_t i=0;i<nRBFs;++i)
				model.rbfs[i].centre = parameter.segment((i+nRBFs+stateDim)*stateDim,stateDim);
			return model;
		}
	};

	template<typename F = RBF<ThinPlateSpline>>
	struct RBFFamily : Family<RBFModel<F>> {
		vector<VectorXd> centres;

		RBFFamily() : RBFFamily(0,0) {}

		RBFFamily( uint32_t dim, uint32_t nRBFs ) : Family<RBFModel<F>>(dim,dim*(dim+nRBFs)), centres(nRBFs) {
			for( auto& c : centres ) {
				c.setZero(dim);
			}
		}

		RBFModel<F> operator()( const VectorXd& parameter ) const {
			RBFModel<F> model(stateDim,(uint32_t)centres.size());
			for(uint32_t i=0;i<stateDim;++i)
				model.linear.col(i) = parameter.segment(i*stateDim,stateDim);
			for(size_t i=0;i<centres.size();++i)
				model.rbfs[i].weight = parameter.segment((i+stateDim)*stateDim,stateDim);
			for(size_t i=0;i<centres.size();++i)
				model.rbfs[i].centre = centres[i];
			return model;
		}

		void ReadText( const char* filename ) {
			ifstream in(filename);
			if( !in )  {
				cout << "RBFFamily::ReadText : file not found " << filename << endl;
				return;
			}
			uint32_t nRBFs;
			in >> stateDim >> nRBFs;
			paramDim = stateDim * ( stateDim + nRBFs );
			centres.resize(nRBFs);
			for( auto& c : centres ) {
				c.setZero(stateDim);
				for(uint32_t j=0;j<stateDim;++j)
					in >> c[j];
			}
		}

		void WriteCSV( const char* filename ) const {
			ofstream out(filename);
			if( !out )  {
				cout << "RBFFamily::WriteCSV : file error " << filename << endl;
				return;
			}
			out << stateDim << "," << centres.size() << endl;

			for( const auto& c : centres ) {
				for(uint32_t j=0;j<stateDim;++j) {
					out << c[j] << ",";
				}
				out << endl;
			}
		}

		MatrixXd ComputeLinear( const VectorXd& x ) const {
			Model model = (*this)( VectorXd::Zero(paramDim) );
			MatrixXd L( stateDim, paramDim );
			L.setZero();
			for(uint32_t i=0;i<stateDim;++i) {
				L.block(0,stateDim*i,stateDim,stateDim).setIdentity() *= x[i];
			}
			for(size_t i=0;i<centres.size();++i) {
				L.block(0,stateDim*(stateDim+i),stateDim,stateDim) = model.rbfs[i].LinearWeight(x);
			}
			return L;
		}

		VectorXd ComputeTranslation( const VectorXd& x ) const {
			VectorXd T(x.size());
			T.setZero();
			return T;
		}

		MatrixXd ComputeLinearDerivative( const VectorXd& x ) const {
			Model model = (*this)( VectorXd::Zero(paramDim) );
			MatrixXd L( stateDim * stateDim, paramDim );
			L.setZero();
			for(uint32_t i=0;i<stateDim;++i) {
				L.block(stateDim*i,stateDim*i,stateDim,stateDim).setIdentity();
			}
			for(size_t i=0;i<centres.size();++i) {
				L.block(0,stateDim*(stateDim+i),stateDim*stateDim,stateDim) = model.rbfs[i].LinearDerivativeWeight(x);
			}
			return L;
		}

		MatrixXd ComputeTranslationDerivative( const VectorXd& x ) const {
			MatrixXd T( x.size(), x.size() );
			T.setZero();
			return T;
		}
	};
}

#endif


