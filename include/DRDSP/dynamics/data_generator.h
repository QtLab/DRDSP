#ifndef INCLUDED_DYNAMICS_DATA_GENERATOR
#define INCLUDED_DYNAMICS_DATA_GENERATOR
#include <iostream>
#include <fstream>
#include <future>
#include "../types.h"
#include "../data/data_system.h"
#include "../dynamics/reduced_data_system.h"

using namespace std;

namespace DRDSP {

	template<typename Family,typename Solver = RKDynamicalSystem<typename Family::Model>>
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
			dtMax(0),
			print(200)
		{}

		explicit DataGenerator( const Family& f ) :
			family(f),
			initial(f.dimension),
			tStart(0),
			tInterval(10),
			dtMax(0),
			print(200)
		{}

		DataGenerator( const Family& f, const State& init ) :
			family(f),
			initial(init),
			tStart(0),
			tInterval(10),
			dtMax(0),
			print(200)
		{}

		void GenerateFiles( const vector<Parameter>& parameters, bool binaryOutput, bool textOutput ) const {
			for( const auto& p : parameters ) {
				GenerateSingleFile(p,binaryOutput,textOutput);
			}
		}

		void GenerateSingleFile( Parameter param, bool binaryOutput, bool textOutput ) const {
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
			Time dtPrint = tInterval / (print-1);
			
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
		}

		DataSet GenerateDataSet( const Parameter& param ) const {
			return GenerateDataSet(param,initial,tInterval);
		}

		DataSet GenerateDataSet( const Parameter& param, Time period ) const {
			return GenerateDataSet(param,initial,period);
		}

		DataSet GenerateDataSet( const Parameter& param, const State& init ) const {
			return GenerateDataSet(param,init,tInterval);
		}

		DataSet GenerateDataSet( const Parameter& param, const State& init, Time period ) const {
			cout << "Parameter " << param << endl;

			Solver solver( family(param) );
			solver.state = init;
			if( dtMax > Time() ) solver.dtMax = dtMax;
			solver.Advance(tStart);

			Time dtPrint = period / (print-1);

			DataSet data( print, family.dimension );

			for(uint32_t i=0;i<print;++i) {
				data[i] = solver.state;
				solver.Advance(dtPrint);
			}
			return data;
		}

		DataSystem GenerateDataSystem( const vector<Parameter>& parameters ) const {
					
			DataSystem data( family.dimension, (uint32_t)parameters.size(), family.parameterDimension );
			
			data.parameters = parameters;

			for(uint32_t i=0;i<parameters.size();++i) {
				data.dataSets[i] = GenerateDataSet( parameters[i] );
			}
			return data;
		}

		DataSystem GenerateDataSystem( const vector<Parameter>& parameters, uint32_t numThreads ) const {
			
			uint32_t numParams = (uint32_t)parameters.size();
			
			DataSystem data( family.dimension, numParams, family.parameterDimension );
			
			data.parameters = parameters;
			vector<future<void>> futures(numThreads);
			
			for(uint32_t i=0;i<numParams;i+=numThreads) {
				uint32_t N = min( numParams - i, numThreads );
				for(uint32_t j=0;j<N;++j) {
					futures[j] = async( launch::async,
						[this]( DataSet& dataSet, const Parameter& parameter ) {
							dataSet = GenerateDataSet( parameter );
						},
						ref(data.dataSets[i+j]), cref(parameters[i+j])
					);
				}
				for(uint32_t j=0;j<N;++j) {
					futures[j].wait();
				}
			}

			return data;
		}

		DataSystem GenerateDataSystem( const vector<Parameter>& parameters, const vector<double>& periods ) {
					
			DataSystem data( family.dimension, (uint32_t)parameters.size(), family.parameterDimension );
			
			data.parameters = parameters;
			
			for(uint16_t i=0;i<parameters.size();++i) {
				data.dataSets[i] = GenerateDataSet( parameters[i], periods[i] );
			}
			return data;
		}

		DataSystem GenerateDataSystem( const vector<Parameter>& parameters, const vector<double>& periods, uint32_t numThreads ) {

			uint32_t numParams = (uint32_t)parameters.size();
			
			DataSystem data( family.dimension, numParams, family.parameterDimension );
			
			data.parameters = parameters;
			vector<future<void>> futures(numThreads);
			
			for(uint32_t i=0;i<numParams;i+=numThreads) {
				uint32_t N = min( numParams - i, numThreads );
				for(uint32_t j=0;j<N;++j) {
					futures[j] = async( launch::async,
						[this]( DataSet& dataSet, const Parameter& parameter, double period ) {
							dataSet = GenerateDataSet( parameter, period );
						},
						ref(data.dataSets[i+j]), cref(parameters[i+j]), periods[i+j]
					);
				}
				for(uint32_t j=0;j<N;++j) {
					futures[j].wait();
				}
			}

			return data;
		}

		DataSystem GenerateUsingInitials( const vector<Parameter>& parameters, const ReducedDataSystem& rdata ) const {
			
			DataSystem data( family.dimension, (uint32_t)parameters.size(), family.parameterDimension );
			
			data.parameters = parameters;
			
			for(uint32_t i=0;i<data.numParameters;++i) {
				data.dataSets[i] = GenerateDataSet( parameters[i], rdata.reducedData[i].points[0] );
			}
			return data;
		}

		DataSystem GenerateUsingInitials( const vector<Parameter>& parameters, const ReducedDataSystem& rdata, uint32_t numThreads ) const {
			
			uint32_t numParams = (uint32_t)parameters.size();
			
			DataSystem data( family.dimension, numParams, family.parameterDimension );
			
			data.parameters = parameters;
			vector<future<void>> futures(numThreads);
			
			for(uint32_t i=0;i<numParams;i+=numThreads) {
				uint32_t N = min( numParams - i, numThreads );
				for(uint32_t j=0;j<N;++j) {
					futures[j] = async( launch::async,
						[this]( DataSet& dataSet, const Parameter& parameter, const VectorXd& init ) {
							dataSet = GenerateDataSet( parameter, init );
						},
						ref(data.dataSets[i+j]), cref(parameters[i+j]), cref(rdata.reducedData[i+j].points[0])
					);
				}
				for(uint32_t j=0;j<N;++j) {
					futures[j].wait();
				}
			}

			return data;
		}

		template<typename F2,typename S2>
		void MatchSettings( const DataGenerator<F2,S2>& gen ) {
			tStart = gen.tStart;
			tInterval = gen.tInterval;
			print = gen.print;
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
