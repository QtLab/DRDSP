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
	struct RBFModel : Model<VectorXd> {
		MatrixXd linear;
		vector<VectorXd> weights;
		vector<RadialFunction<F>> rbfs;
		uint32_t numRBFs;
		mt19937 mt;
		uniform_real_distribution<double> dist;

		RBFModel() : numRBFs(0) {}
		
		RBFModel( uint32_t dim, uint32_t nRBFs ) :
			Model<VectorXd>(dim),
			numRBFs(nRBFs),
			weights(nRBFs),
			rbfs(nRBFs)
		{
			linear.setZero(dimension,dimension);
			for(uint32_t i=0;i<numRBFs;++i) {
				weights[i].setZero(dimension);
				rbfs[i].centre.setZero(dimension);
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

		void SetCentresRandom( const AABB& box ) {
			VectorXd diff = box.bMax - box.bMin;
			for(uint32_t i=0;i<numRBFs;++i)
				for(uint32_t j=0;j<dimension;++j) {
					rbfs[i].centre(j) = box.bMin(j) + diff(j) * dist(mt);
				}
		}

		void LoadCentresText( const char* filename ) {
			ifstream in(filename);
			if( !in ) return;

			for(uint32_t k=0;k<numRBFs;++k)
				for(uint32_t j=0;j<dimension;++j)
					in >> rbfs[k].centre(j);
		}

		void LoadCentresBinary( const char* filename ) {
			ifstream in(filename);
			if( !in ) return;

			for(uint32_t k=0;k<numRBFs;++k)
				in.read( (char*)&rbfs[k].centre(0), sizeof(double)*dimension );
		}

		void WriteCSV( const char *filename ) const {
			ofstream out(filename);
			out.precision(16);
			out << dimension << "," << numRBFs << endl;
			for(uint32_t i=0;i<dimension;++i) {	
				for(uint32_t j=0;j<dimension;++j)
					out << linear(i,j) << ",";
				out << endl;
			}
			out << endl;
			for(uint32_t i=0;i<dimension;++i) {	
				for(uint32_t j=0;j<numRBFs;++j)
					out << weights[j](i) << ",";
				out << endl;
			}
			out << endl;
			for(uint32_t i=0;i<dimension;++i) {	
				for(uint32_t j=0;j<numRBFs;++j)
					out << rbfs[j].centre(i) << ",";
				out << endl;
			}
		}

	};

}

#endif
