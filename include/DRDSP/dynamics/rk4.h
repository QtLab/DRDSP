#ifndef INCLUDED_DYNAMICS_RK4
#define INCLUDED_DYNAMICS_RK4

namespace DRDSP {
	template<typename F>
	struct RK4 {
		typedef double Time;
		typedef typename F::State State;
		F f;

		explicit RK4( const F& f ) : f(f) {}

		void Integrate( State& x, double& t, double dt ) {
			double halfdt = 0.5 * dt;

			RKt[0] = f( x, t );
			RKt[1] = f( x + halfdt * RKt[0], t + halfdt );
			RKt[2] = f( x + halfdt * RKt[1], t + halfdt );
			RKt[3] = f( x + dt * RKt[2], t + dt );

			x += ( RKt[0] + 2.0 * RKt[1] + 2.0 * RKt[2] + RKt[3] ) * (dt / 6.0);
			t += dt;
		}

	protected:
		State RKt[4];
	};
}

#endif
