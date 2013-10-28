#ifndef INCLUDED_DYNAMICS_DYNAMICALSYSTEM
#define INCLUDED_DYNAMICS_DYNAMICALSYSTEM
#include "../types.h"
#include "rk.h"

namespace DRDSP {
	template<typename TTime,typename TState>
	struct DynamicalSystem {
		TState state;
		virtual void Advance( TTime dt ) = 0;
	};

	template<typename TTime,typename TState>
	struct ContinuousDynamicalSystem : DynamicalSystem<TTime,TState> {
		Solver<TTime,TState>& solver;
		double dtmax;
		
		ContinuousDynamicalSystem( Solver<TTime,TState>& S ) : solver(S), dtmax(0) {}
		
		void Advance( TTime dt ) {
			double t = 0.0, dta = dt;
			uint32_t n = 1;
			if( dt > dtmax && dtmax > 0.0 ) {
				n = uint32_t(dt / dtmax + 1.0);
				dta = dt / n;
			}
			for(uint32_t i=0;i<n;i++)
				solver.Integrate(state,t,dta);
		}

	};

	template<typename TTime,typename TState>
	struct RKDynamicalSystem : ContinuousDynamicalSystem<TTime,TState> {
		RKDynamicalSystem( SolverFunction<TTime,TState>& F ) : rk(F), ContinuousDynamicalSystem<TTime,TState>(rk) {}
	protected:
		RK<TState> rk;
	}; 

}

#endif
