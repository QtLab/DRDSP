#ifndef INCLUDED_DYNAMICS_RBF_MODEL
#define INCLUDED_DYNAMICS_RBF_MODEL
#include "../types.h"
#include "model.h"
#include "radial_basis.h"

using namespace std;

namespace DRDSP {

	template<typename F = RBF<ThinPlateSpline>>
	struct RBFModel : Model<> {
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
