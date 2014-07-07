#ifndef INCLUDED_DYNAMICS_RBF_MODEL
#define INCLUDED_DYNAMICS_RBF_MODEL
#include "../types.h"
#include "model.h"
#include "radial_basis.h"
#include "../data/aabb.h"
#include <random>

using namespace std;

namespace DRDSP {

	template<typename F = ThinPlateSpline>
	struct RBFModel2 : Model<> {
		MatrixXd linear;
		vector<VectorXd> weights;
		vector<RBF<F>> rbfs;
		uint32_t numRBFs;

		RBFModel2() : numRBFs(0) {}
		
		RBFModel2( uint32_t dim, uint32_t nRBFs ) :
			Model<>(dim),
			numRBFs(nRBFs),
			weights(nRBFs),
			rbfs(nRBFs)
		{
			linear.setZero(stateDim,stateDim);
			for(uint32_t i=0;i<numRBFs;++i) {
				weights[i].setZero(stateDim);
				rbfs[i].centre.setZero(stateDim);
			}
		}

		VectorXd operator()( const VectorXd& x ) const {
			VectorXd sum = linear * x;
			for(uint32_t i=0;i<numRBFs;++i)
				sum += weights[i] * rbfs[i](x);
			return sum;
		}

		MatrixXd Partials( const VectorXd& x ) const {
			MatrixXd sum = linear;
			for(uint32_t i=0;i<numRBFs;++i)
				sum += weights[i] * rbfs[i].Derivative(x).transpose();
			return sum;
		}

		void LoadCentresText( const char* filename ) {
			ifstream in(filename);
			if( !in ) return;

			for(uint32_t k=0;k<numRBFs;++k)
				for(uint32_t j=0;j<stateDim;++j)
					in >> rbfs[k].centre(j);
		}

		void LoadCentresBinary( const char* filename ) {
			ifstream in(filename);
			if( !in ) return;

			for(uint32_t k=0;k<numRBFs;++k)
				in.read( (char*)&rbfs[k].centre(0), sizeof(double)*stateDim );
		}

		void WriteCSV( const char *filename ) const {
			ofstream out(filename);
			out.precision(16);
			out << dimension << "," << numRBFs << endl;
			for(uint32_t i=0;i<stateDim;++i) {	
				for(uint32_t j=0;j<stateDim;++j)
					out << linear(i,j) << ",";
				out << endl;
			}
			out << endl;
			for(uint32_t i=0;i<stateDim;++i) {	
				for(uint32_t j=0;j<numRBFs;++j)
					out << weights[j](i) << ",";
				out << endl;
			}
			out << endl;
			for(uint32_t i=0;i<stateDim;++i) {	
				for(uint32_t j=0;j<numRBFs;++j)
					out << rbfs[j].centre(i) << ",";
				out << endl;
			}
		}

	};

	template<typename F = RBF<ThinPlateSpline>>
	struct RBFModel : Model<> {
		typedef typename F::RadialType RadialType;
		MatrixXd linear;
		vector<F> rbfs;

		RBFModel() = default;
		
		RBFModel( uint32_t dim, uint32_t nRBFs ) :
			Model<>(dim),
			rbfs(nRBFs)
		{
			linear.setZero(stateDim,stateDim);
			for( auto& r : rbfs ) {
				r.weight.setZero(stateDim);
				r.centre.setZero(stateDim);
			}
		}

		VectorXd operator()( const VectorXd& x ) const {
			VectorXd sum = linear * x;
			for( const auto& r : rbfs )
				sum += r(x);
			return sum;
		}

		MatrixXd Partials( const VectorXd& x ) const {
			MatrixXd sum = linear;
			for( const auto& r : rbfs )
				sum += r.Derivative(x);
			return sum;
		}

		void LoadCentresText( const char* filename ) {
			ifstream in(filename);
			if( !in ) return;

			for( auto& r : rbfs )
				for(uint32_t j=0;j<stateDim;++j)
					in >> r.centre[j];
		}

		void LoadCentresBinary( const char* filename ) {
			ifstream in(filename);
			if( !in ) return;

			for( auto& r : rbfs )
				in.read( (char*)&r.centre[0], sizeof(double)*stateDim );
		}

		void WriteCSV( const char* filename ) const {
			ofstream out(filename);
			out.precision(16);
			out << dimension << "," << numRBFs << endl;
			for(uint32_t i=0;i<stateDim;++i) {	
				for(uint32_t j=0;j<stateDim;++j)
					out << linear(i,j) << ",";
				out << endl;
			}
			out << endl;
			for( const auto& r : rbfs ) {
				for(uint32_t i=0;i<stateDim;++i) {
					out << r.weight[i] << ",";
				}
				out << endl;
			}
			out << endl;
			for( const auto& r : rbfs ) {
				for(uint32_t i=0;i<stateDim;++i) {
					out << r.centre[i] << ",";
				}
				out << endl;
			}
		}

	};
}

#endif
