#ifndef INCLUDED_DYNAMICS_RBF_FAMILY
#define INCLUDED_DYNAMICS_RBF_FAMILY
#include <iostream>
#include <fstream>
#include "../types.h"
#include "affineParameterMap.h"
#include "rbf_model.h"
#include <Eigen/LU>

using namespace std;

namespace DRDSP {
	template<typename F = ThinPlateSpline>
	struct RBFFamily : Family<RBFModel<F>,VectorXd> {
		AffineParameterMap affine;
		RBFModel<F> model;

		RBFFamily() = default;

		RBFFamily( uint32_t dim, uint32_t paramDim, uint32_t nRBFs ) : 
			Family<RBFModel<F>,VectorXd>( dim, paramDim ),
			model(dim,nRBFs),
			affine(dimension,nRBFs,paramDim)
		{}

		RBFModel<F> operator()( const VectorXd& parameter ) const {
			MatrixXd z = affine(parameter);

			RBFModel<F> r = model;

			r.linear = z.block(0,0,dimension,dimension);

			for(uint32_t i=0;i<r.numRBFs;i++) {
				r.weights[i] = z.col(dimension+i);
			}
			return r;
		}

		VectorXd VectorField( const VectorXd& x, const VectorXd &parameter ) {
			MatrixXd z = affine(parameter);

			model.linear = z.block(0,0,dimension,dimension);

			for(uint32_t i=0;i<model.numRBFs;i++) {
				model.weights[i] = z.col(dimension+i);
			}
			return model(x);
		}

		MatrixXd Partials( const VectorXd& x, const VectorXd& parameter ) {
			MatrixXd z = affine(parameter);

			model.linear = z.block(0,0,dimension,dimension);

			for(uint32_t i=0;i<model.numRBFs;i++) {
				model.weights[i] = z.col(dimension+i);
			}
			return model.Partials(x);
		}

		void WriteCSV( const char *filename ) const {
			ofstream out(filename);
			out.precision(16);
			out << dimension << "," << model.numRBFs << "," << parameterDimension << endl;
			for(int i=0;i<affine.coeffs.rows();i++) {	
				for(int j=0;j<affine.coeffs.cols();j++)
					out << affine.coeffs(i,j) << ",";
				out << endl;
			}
			for(uint32_t k=0;k<model.numRBFs;k++) {
				for(uint32_t j=0;j<dimension;j++)
					out << model.rbfs[k].centre(j) << ",";
				out << endl;
			}
		}

		void ReadText( const char* filename ) {
			ifstream in(filename);
			if( !in )  {
				cout << "RBFFamily::ReadText : file not found " << filename << endl;
				return;
			}
			uint32_t numRBFs;
			
			in >> dimension >> numRBFs >> parameterDimension;
			
			affine = AffineParameterMap( dimension, numRBFs, parameterDimension );
			model = RBFModel<F>( dimension, numRBFs );

			for(uint32_t i=0;i<affine.coeffs.rows();i++)
				for(uint32_t j=0;j<affine.coeffs.cols();j++)
					in >> affine.coeffs(i,j);
			for(uint32_t i=0;i<numRBFs;i++)
				for(uint32_t j=0;j<dimension;j++)
					in >> model.rbfs[i].centre(j);
		}

	};
}

#endif




