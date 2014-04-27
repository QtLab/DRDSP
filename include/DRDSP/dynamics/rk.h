#ifndef INCLUDED_DYNAMICS_RK
#define INCLUDED_DYNAMICS_RK
#include <stdint.h>

namespace DRDSP {
	const double RK1_Table[] = {1.0,0.0};
	const double RK2_Table[] = {0.5,0.5,1.0,1.0,0.0};
	const double RK3_Table[] = {1.0/6.0,1.0/6.0,1.0/6.0,1.0,-1.0,2.0,0.5,0.5,0.0};
	const double RK4_Table[] = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0,1.0,0.0,0.0,1.0,0.5,0.0,0.5,0.5,0.5,0.0};

	template<typename F,uint8_t order=4>
	struct RK {
		typedef double Time;
		typedef typename F::State State;
		F f;

		explicit RK( const F& f ) : f(f) {
			
			for(uint8_t i=0;i<order-1;++i)
				for(uint8_t j=0;j<order-1;++j)
					RK_a[i][j] = 0.0;
			
			LoadButcherTableau();
		}

		void Integrate( State& x, double& t, double dt ) {

			RKt[0] = f(x,t);

			double u2;
			for(uint16_t n=1;n<order;++n) {
				u2 = dt * RK_a[n-1][0];
				RKTemp = x + RKt[0] * u2;
				
				for(uint16_t m=1;m<n;++m) {
					u2 = dt * RK_a[n-1][m-1];
					RKTemp += RKt[m] * u2;
				}
				u2 = t + dt * RK_c[n];
				RKt[n] = f(RKTemp,u2);
			}

			RKTemp = RKt[0] * RK_b[0];
			for(uint8_t n=1;n<order;++n)
				RKTemp += RKt[n] * RK_b[n];

			x += RKTemp * dt;
			t += dt;
		}

	protected:
		double RK_a[order-1][order-1];
		double RK_b[order];
		double RK_c[order];
		State RKt[order], RKTemp;

		bool LoadButcherTableau() {
			const double *table;
			switch(order) {
				case 1:
					table = RK1_Table;
				break;
				case 2:
					table = RK2_Table;
				break;
				case 3:
					table = RK3_Table;
				break;
				case 4:
					table = RK4_Table;
				break;
				default:
				return false;
			};

			uint16_t k=0;
			for(uint8_t i=0;i<order;++i) {
				RK_b[i] = table[k++];
			}

			for(uint8_t j=0;j<order;++j) {
				RK_c[j] = table[k++];
				for(uint8_t i=0;i<order-1-j;++i) {
					RK_a[order-2-j][i] = table[k++];
				}
			}
			return true;
		}

	};

}

#endif


