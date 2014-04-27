#ifndef INCLUDED_DYNAMICS_GENERATE_DATA
#define INCLUDED_DYNAMICS_GENERATE_DATA
#include <iostream>
#include <fstream>
#include "../types.h"
#include "../data/data_system.h"
#include "../dynamics/reduced_data_system.h"

using namespace std;

namespace DRDSP {

	template<typename Family,typename Solver = RKDynamicalSystem<SolverFunctionFromModel<typename Family::Model>>>
	struct DataGenerator {
		typedef typename Family::Parameter Parameter;
		typedef typename Family::Model Model;
		typedef typename Model::State State;
		typedef typename Solver::Time Time;
		Family family;
		State initial;
		Time tStart, tInterval, dtMax;
		uint32_t print;
		bool binaryOutput, textOutput;
	
		DataGenerator() :
			tStart(0),
			tInterval(10),
			dtMax(0.001),
			print(200),
			binaryOutput(true),
			textOutput(false)
		{}

		explicit DataGenerator( const Family& f ) :
			family(f),
			initial(f.dimension),
			tStart(0),
			tInterval(10),
			dtMax(0.001),
			print(200),
			binaryOutput(true),
			textOutput(false)
		{}

		DataGenerator( const Family& f, const State& init ) :
			family(f),
			initial(init),
			tStart(0),
			tInterval(10),
			dtMax(0.001),
			print(200),
			binaryOutput(true),
			textOutput(false)
		{}

		void GenerateFiles( const vector<Parameter>& parameters ) const {
			for( const auto& p : parameters ) {
				GenerateSingleFile(p);
			}
		}

		void GenerateSingleFile( Parameter param ) const {
			ofstream outBin, outTxt;
			stringstream fn;
			outTxt.precision(16);

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

			Solver solver( SolverFunctionFromModel<Model>( family(param) ) );
			solver.state = initial;
			solver.dtMax = dtMax;
			solver.Advance(tStart);

			Time t = 0.0;
			Time dtPrint = tInterval / print;
			
			if( binaryOutput ) {
				outBin.write((char*)&print,sizeof(uint32_t));
				outBin.write((char*)&param,sizeof(Parameter));
			}
			if( textOutput ) outTxt << print << "," << param << endl;

			while( t <= tInterval ) {
			
				if( binaryOutput ) {
					outBin.write((const char*)&solver.state(0),sizeof(double)*solver.state.size());
				}
				if( textOutput ) {
					for(int i=0;i<solver.state.size();++i) {
						outTxt << solver.state(i) << ",";
					}
					outTxt << endl;
				}
				solver.Advance(dtPrint);
				t += dtPrint;
			}
			if( binaryOutput ) outBin.close();
			if( textOutput ) outTxt.close();
		}

		DataSet GenerateDataSet( const Parameter& param ) const {
			return GenerateDataSet(param,initial);
		}
		
		DataSet GenerateDataSet( const Parameter& param, const State& init ) const {
			cout << "Parameter " << param << endl;

			Solver solver( SolverFunctionFromModel<Model>( family(param) ) );
			solver.state = init;
			solver.dtMax = dtMax;
			solver.Advance(tStart);

			Time dtPrint = tInterval / print;

			DataSet data( print, family.dimension );

			for(uint32_t i=0;i<print;i++) {
				data[i] = solver.state;
				solver.Advance(dtPrint);
			}
			return data;
		}

		DataSystem GenerateDataSystem( const vector<Parameter>& parameters ) const {
					
			DataSystem data( family.dimension, (uint32_t)parameters.size(), family.parameterDimension );
			
			for(uint16_t i=0;i<parameters.size();i++) {
				data.parameters[i] = parameters[i];
				data.dataSets[i] = GenerateDataSet( parameters[i] );
			}

			return data;
		}

		DataSystem GenerateUsingInitials( const vector<Parameter>& parameters, const ReducedDataSystem& rdata ) const {
			
			DataSystem data( family.dimension, (uint32_t)parameters.size(), family.parameterDimension );
			
			data.parameters = parameters;
			
			for(uint16_t i=0;i<parameters.size();i++) {
				data.dataSets[i] = GenerateDataSet( parameters[i], rdata.reducedData[i].points[0] );
			}

			return data;
		}

		template<typename F2,typename S2>
		void MatchSettings( const DataGenerator<F2,S2>& gen ) {
			tStart = gen.tStart;
			tInterval = gen.tInterval;
			print = gen.print;
			binaryOutput = gen.binaryOutput;
			textOutput = gen.textOutput;
			initial = gen.initial;
		}

	};

	vector<VectorXd> ParameterList( const VectorXd& pMin, const VectorXd& pMax, uint32_t N ) {
		vector<VectorXd> v(N);
		VectorXd pDelta = ( pMax - pMin ) / ( N - 1 );
		v[0] = pMin;
		for(uint32_t i=1;i<N;++i) {
			v[i] = v[i-1] + pDelta;
		}
		return v;
	}

	vector<VectorXd> ParameterList( double pMin, double pMax, uint32_t N ) {
		vector<VectorXd> v(N);
		double pDelta = ( pMax - pMin ) / ( N - 1 );
		VectorXd temp(1);
		temp[0] = pMin;
		v[0] = temp;
		for(uint32_t i=1;i<N;++i) {
			temp[0] = pDelta;
			v[i] = v[i-1] + temp;
		}
		return v;
	}

}

#endif
