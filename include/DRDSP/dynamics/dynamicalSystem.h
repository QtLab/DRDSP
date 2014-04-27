#ifndef INCLUDED_DYNAMICS_DYNAMICALSYSTEM
#define INCLUDED_DYNAMICS_DYNAMICALSYSTEM
#include <stdint.h>
#include "rk.h"

namespace DRDSP {

	template<typename State,typename Map>
	struct DiscreteDynamicalSystem {
		State state;
		Map map;

		DiscreteDynamicalSystem() = default;
		
		explicit DiscreteDynamicalSystem( const Map& map ) : map(map) {}

		explicit DiscreteDynamicalSystem( const State& state ) : state(state) {}

		DiscreteDynamicalSystem( const Map& map, const State& state ) : map(map), state(state) {}

		void Advance( uint32_t dt = 1 ) {
			for(uint32_t i=0;i<dt;++i)
				state = map(state);
		}

	};

	template<typename Solver>
	struct ContinuousDynamicalSystem {
		typedef double Time;
		typedef typename Solver::State State;
		State state;
		Solver solver;
		double dtMax;
		
		explicit ContinuousDynamicalSystem( const Solver& solver ) : solver(solver), dtMax(0) {}

		ContinuousDynamicalSystem( const Solver& solver, const State& state ) : solver(solver), state(state), dtMax(0) {}

		void Advance( double dt ) {
			double t = 0.0, dta = dt;
			uint32_t n = 1;
			if( dt > dtMax && dtMax > 0.0 ) {
				n = uint32_t(dt / dtMax + 1.0);
				dta = dt / n;
			}
			for(uint32_t i=0;i<n;++i)
				solver.Integrate(state,t,dta);
		}

	};

	template<typename F>
	struct RKDynamicalSystem : ContinuousDynamicalSystem<RK<F>> {
		explicit RKDynamicalSystem( const F& f ) : ContinuousDynamicalSystem<RK<F>>(RK<F>(f)) {}
	};

}

#endif



