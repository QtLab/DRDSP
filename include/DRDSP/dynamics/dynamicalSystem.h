#ifndef INCLUDED_DYNAMICS_DYNAMICALSYSTEM
#define INCLUDED_DYNAMICS_DYNAMICALSYSTEM
#include <stdint.h>
#include "rk4.h"

namespace DRDSP {

	template<typename State,typename Map>
	struct DiscreteDynamicalSystem {
		State state;
		Map map;
		
		DiscreteDynamicalSystem() = default;
		
		explicit DiscreteDynamicalSystem( const Map& map ) : map(map) {}

		DiscreteDynamicalSystem( const Map& map, const State& state ) : map(map), state(state) {}

		void Advance( uint32_t dt = 1 ) {
			for(uint32_t i=0;i<dt;++i)
				state = map(state);
		}
	};

	struct IdentityWrapFunction {
		template<typename T>
		void operator()( T& ) const {}
	};

	template<typename Solver,typename WrapFunction = IdentityWrapFunction>
	struct ContinuousDynamicalSystem {
		typedef double Time;
		typedef typename Solver::State State;
		State state;
		Solver solver;
		WrapFunction wrap;
		double dtMax = 0.0;
		
		explicit ContinuousDynamicalSystem( const Solver& solver ) : solver(solver) {}

		ContinuousDynamicalSystem( const Solver& solver, const State& state ) : solver(solver), state(state) {}

		void Advance( double dt ) {
			double dta = dt;
			uint32_t n = 1;
			if( dt > dtMax && dtMax > 0.0 ) {
				n = uint32_t(dt / dtMax + 1.0);
				dta = dt / n;
			}
			for(uint32_t i=0;i<n;++i) {
				solver.Integrate(state,dta);
				wrap(state);
			}
		}
	};

	template<typename F,typename WrapFunction = IdentityWrapFunction>
	struct RKDynamicalSystem : ContinuousDynamicalSystem<RK4<F>,WrapFunction> {
		explicit RKDynamicalSystem( const F& f ) : ContinuousDynamicalSystem<RK4<F>,WrapFunction>(RK4<F>(f)) {}
	};

}

#endif

