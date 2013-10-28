#ifndef INCLUDED_DYNAMICS_GENERATE_DATA
#define INCLUDED_DYNAMICS_GENERATE_DATA
#include <iostream>
#include <fstream>
#include "../types.h"
#include "dynamicalSystem.h"
#include "model_orig.h"

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

	struct GenerateData {

		GenerateData( ModelParameterized& m ) : model(m), rk(model), tStart(0), tEnd(10), print(200), pMin(0), pMax(1), pDelta(0.1), binaryOutput(true), textOutput(false) {
			initial.setZero(m.dimension);
		}

		double pMin, pMax, pDelta;
		double tStart, tEnd;
		uint32_t print;
		bool binaryOutput, textOutput;
		VectorXd initial;

		ModelParameterizedInterface model;
		RKDynamicalSystem<double,VectorXd> rk;

		void Generate() {
			for(double p=pMin;p<=pMax;p+=pDelta) {
				GenerateSingle(p);
			}
		}

		void GenerateSingle( double param ) {
			ofstream outBin, outTxt;
			stringstream fn;
			outTxt.precision(16);
			model.parameter(0) = param;
			fn.str("");
			fn << "data/pp" << param;
			cout << fn.str() << endl;
		
			if( binaryOutput ) {
				fn << ".bin";
				outBin.open(fn.str(),ios::binary);
			}
			if( textOutput ) {
				fn.str("");
				fn << "data/pp" << param << ".csv";
				outTxt.open(fn.str());
			}

			rk.state = initial;
			rk.Advance(tStart);

			double t = tStart;
			double dtPrint = (tEnd-tStart) / print;
			
			if( binaryOutput ) {
				outBin.write((char*)&print,sizeof(uint32_t));
				outBin.write((char*)&param,sizeof(double));
			}
			if( textOutput )
				outTxt << print << "," << param << endl;

			while( t <= tEnd ) {
			
				if( binaryOutput ) {
					outBin.write((const char*)&rk.state(0),sizeof(double)*rk.state.size());
				}
				if( textOutput ) {
					for(int i=0;i<rk.state.size();i++) {
						outTxt << rk.state(i) << ",";
					}
					outTxt << endl;
				}
				rk.Advance(dtPrint);
				t += dtPrint;
			}
			if( binaryOutput ) outBin.close();
			if( textOutput ) {
				outTxt.close();
			}
		}

	};

}

#endif

