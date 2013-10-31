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

	template<typename TState>
	struct WrapFunction {
		virtual void operator()( TState& x ) const {}
		static const WrapFunction<TState> identity;
	};

	template<typename TTime,typename TState>
	struct ContinuousDynamicalSystem : DynamicalSystem<TTime,TState> {
		Solver<TTime,TState>& solver;
		const WrapFunction<TState>& wrap;
		double dtMax;
		
		explicit ContinuousDynamicalSystem( Solver<TTime,TState>& S ) : solver(S), wrap(WrapFunction<TState>::identity), dtMax(0) {}

		ContinuousDynamicalSystem( Solver<TTime,TState>& S, const WrapFunction<TState>& W ) : solver(S), wrap(W), dtMax(0) {}

		void Advance( TTime dt ) {
			double t = 0.0, timeStep = dt;
			uint32_t n = 1;
			if( dt > dtMax && dtMax > 0.0 ) {
				n = uint32_t(dt / dtMax + 1.0);
				timeStep = dt / n;
			}
			for(uint32_t i=0;i<n;i++) {
				solver.Integrate(state,t,timeStep);
				wrap(state);
			}
		}

	};

	template<typename TTime,typename TState>
	struct RKDynamicalSystem : ContinuousDynamicalSystem<TTime,TState> {
		explicit RKDynamicalSystem( SolverFunction<TTime,TState>& F ) : rk(F), ContinuousDynamicalSystem<TTime,TState>(rk) {}
		RKDynamicalSystem( SolverFunction<TTime,TState>& F, const WrapFunction<TState>& W ) : rk(F), ContinuousDynamicalSystem<TTime,TState>(rk,W) {}
	protected:
		RK<TState> rk;
	}; 

}

#endif
