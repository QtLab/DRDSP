#ifndef INCLUDED_DYNAMICS_BIFURCATION
#define INCLUDED_DYNAMICS_BIFURCATION
#include <fstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <stdint.h>
#include "../bitmap.h"

using namespace std;

namespace DRDSP {

	template<typename Parameter,typename Value>
	struct BifurcationDiagram {
		vector<pair<Parameter,Value>> data;

		BifurcationDiagram() = default;

		explicit BifurcationDiagram( uint32_t numParms ) {
			data.reserve( numParms );
		}

		void WriteCSV( const char* filename ) const {
			ofstream out(filename);
			Parameter param = data[0].first;
			out << param << ",";
			for(const auto& p : data) {
				if( p.first != param ) {
					out << endl << p.first << ",";
					param = p.first;
				}
				out << p.second << ",";
			}
		}

		void WriteBitmap( const char* filename, uint32_t Ny ) const {			
			
			auto minmax = minmax_element( begin(data), end(data),
				[](const pair<Parameter,Value>& x, const pair<Parameter,Value>& y){ return x.second < y.second; }
			);

			Value range = minmax.second->second - minmax.first->second;

			Bitmap bmp( CountParameters(), Ny );
			bmp.Clear(255,255,255);

			uint32_t i = 0;
			Parameter param = data[0].first;

			for(const auto& p : data) {
				if( p.first != param ) {
					++i;
					param = p.first;
				}
				uint32_t k = Ny-1 - uint32_t(( p.second - minmax.first->second )/range * (Ny-1));
				bmp(i,k).r = 0;
				bmp(i,k).g = 0;
				bmp(i,k).b = 255;
			}

			bmp.WriteFile(filename);
		}

		uint32_t CountParameters() const {
			if( data.size() == 0 ) return 0;
			uint32_t i = 1;
			Parameter param = data[0].first;

			for(const auto& p : data) {
				if( p.first != param ) {
					++i;
					param = p.first;
				}
			}
			return i;
		}

	};

	template<typename Family>
	struct BifurcationDiagramGenerator {
		typedef double Time;
		typedef typename Family::Model Model;
		typedef typename Model::State State;
		typedef typename Family::Parameter Parameter;
		typedef double Value;

		uint32_t numParams;
		Parameter minParam, maxParam;
		Time startTime, endTime, dt;
		State initial;
		
		template<typename Condition,typename GetValue>
		BifurcationDiagram<Parameter,Value> Generate( Family& family, Condition&& condition, GetValue&& getValue ) const {
			State prevState, newState;

			Parameter dParam = ( maxParam - minParam ) / ( numParams - 1 );

			BifurcationDiagram<Parameter,Value> bifurcationDiagram( numParams );

			Parameter parameter = minParam;
			for(uint32_t i=0;i<numParams;++i) {
				RKDynamicalSystem<SolverFunctionFromModel<Model>> system( SolverFunctionFromModel<Model>( family(parameter) ) );
				system.state = initial;
				system.dtMax = 0.0005;
				system.Advance(startTime);
				prevState = system.state;
				for(Time t=startTime;t<endTime;t+=dt) {
					system.Advance(dt);
					newState = system.state;
					if( condition(prevState,newState) ) {
						bifurcationDiagram.data.emplace_back( parameter, getValue(prevState,newState) );
					}
					prevState = newState;
				}
				parameter += dParam;
			}
			return bifurcationDiagram;
		}


	};

}

#endif

