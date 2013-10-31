#ifndef INCLUDED_DYNAMICS_RK
#define INCLUDED_DYNAMICS_RK
#include "../types.h"
#include "solver.h"
#include "../misc.h"

namespace DRDSP {
	const double RK1_Table[] = {1.0,0.0};
	const double RK2_Table[] = {0.5,0.5,1.0,1.0,0.0};
	const double RK3_Table[] = {1.0/6.0,1.0/6.0,1.0/6.0,1.0,-1.0,2.0,0.5,0.5,0.0};
	const double RK4_Table[] = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0,1.0,0.0,0.0,1.0,0.5,0.0,0.5,0.5,0.5,0.0};

	template<typename TState>
	struct RK : Solver<double,TState> {

		explicit RK( SolverFunction<double,TState>& f ) : Solver(f), created(false) { Create(4); }

		RK( uint8_t RK_order, const SolverFunction<double,TState>& f ) : Solver(f), created(false), order(4) {
			Create(RK_order);
		}

		~RK() { Destroy(); }

		void Create( uint8_t RK_order ) {
			if( created ) Destroy();
			order = Clamp<uint8_t>(RK_order,1,4);
			RKt = new TState [order];
			RK_c = new double [order];
			RK_b = new double [order];
			RK_a = new double *[order-1];
			for(uint8_t i=0;i<order-1;i++)
				RK_a[i] = new double [order-1];

			created = true;

			for(uint8_t i=0;i<order-1;i++)
				for(uint8_t j=0;j<order-1;j++)
					RK_a[i][j] = 0.0;

			LoadButcherTableau();
		}

		void Destroy() {
			if( !created ) return;
			delete[] RKt;
			for(uint8_t i=0;i<order-1;i++)
				delete[] RK_a[i];
			delete[] RK_a;
			delete[] RK_b;
			delete[] RK_c;
			created = false;
		}

		void Integrate( TState& x, double& t, double dt ) {

			RKt[0] = Function(x,t);

			double u2;
			for(uint16_t n=1;n<order;n++) {
				u2 = dt * RK_a[n-1][0];
				RKTemp = x + RKt[0] * u2;
				
				for(uint16_t m=1;m<n;m++) {
					u2 = dt * RK_a[n-1][m-1];
					RKTemp += RKt[m] * u2;
				}
				u2 = t + dt * RK_c[n];
				RKt[n] = Function(RKTemp,u2);
			}

			RKSum = RKt[0] * RK_b[0];
			for(uint8_t n=1;n<order;n++)
				RKSum += RKt[n] * RK_b[n];

			x += RKSum * dt;
			t += dt;
		}

	protected:
		bool created;
		uint8_t order;
		double **RK_a, *RK_b, *RK_c;
		TState *RKt, RKTemp, RKSum;

		void LoadButcherTableau() {
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
				default:
				case 4:
					table = RK4_Table;
				break;
			};

			uint16_t k=0;
			for(uint8_t i=0;i<order;i++) {
				RK_b[i] = table[k++];
			}

			for(uint8_t j=0;j<order;j++) {
				RK_c[j] = table[k++];
				for(uint8_t i=0;i<order-1-j;i++) {
					RK_a[order-2-j][i] = table[k++];
				}
			}
		}

	};

}

#endif


