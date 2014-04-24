#ifndef INCLUDED_DYNAMICS_GENERATE_DATA
#define INCLUDED_DYNAMICS_GENERATE_DATA
#include <iostream>
#include <fstream>
#include "../types.h"
#include "dynamicalSystem.h"
#include "model.h"
#include "../data/data_set.h"

using namespace std;

namespace DRDSP {

	struct ModelInterface : SolverFunction<double,VectorXd> {
		Model& model;

		explicit ModelInterface( Model& m ) : model(m) {}

		VectorXd operator()( const VectorXd& x, double t ) {
			return model.VectorField(x);
		}
	};

	struct ModelParameterizedInterface : SolverFunction<double,VectorXd> {
		ModelParameterized& model;
		VectorXd parameter;

		explicit ModelParameterizedInterface( ModelParameterized& m ) : model(m) {
			parameter.setZero(model.parameterDimension);
		}

		VectorXd operator()( const VectorXd& x, double t ) {
			return model.VectorField(x,parameter);
		}
	};

	struct DataGenerator {

		double pMin, pMax, pDelta;
		double tStart, tInterval;
		uint32_t print;
		bool binaryOutput, textOutput;
		VectorXd initial;
		ModelParameterizedInterface model;
		RKDynamicalSystem<double,VectorXd> rk;

		DataGenerator( ModelParameterized& m );
		void GenerateSingleFile( double param );
		void GenerateFiles();
		DataSet GenerateDataSet( double param );
		DataSystem GenerateDataSystem();
		void MatchSettings( const DataGenerator& gen );
	};





}

#endif

