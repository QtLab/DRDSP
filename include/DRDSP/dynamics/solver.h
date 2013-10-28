#ifndef INCLUDED_DYNAMICS_SOLVER
#define INCLUDED_DYNAMICS_SOLVER

namespace DRDSP {

	template<typename TTime,typename TState>
	struct SolverFunction {
		virtual TState operator()( const TState& x, TTime t ) = 0;
	};

	template<typename TTime,typename TState>
	struct Solver {
		SolverFunction<TTime,TState>& Function;

		Solver( SolverFunction<TTime,TState>& f ) : Function(f) {}
		virtual void Integrate( TState& x, TTime& t, TTime dt ) = 0;
	};

	template<typename TTime,typename TState>
	struct Euler : Solver<TTime,TState> {
		void Integrate( TState& x, TTime& t, TTime dt ) {
			x += Function(x,t) * dt;
			t += dt;
		}
	};

}

#endif


