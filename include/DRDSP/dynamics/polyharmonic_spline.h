#ifndef INCLUDED_DYNAMICS_POLYHARMONIC_SPLINE
#define INCLUDED_DYNAMICS_POLYHARMONIC_SPLINE
#include <cmath>

namespace DRDSP {
	
	using std::log;

	template<int N>
	struct PolyharmonicSpline {
		template<typename T>
		T operator()( T r ) const {
			return r * r * PolyharmonicSpline<N-2>()(r);
		}
		template<typename T>
		T Derivative( T r ) const {
			return r * ( T(2) * PolyharmonicSpline<N-2>()(r) + r * PolyharmonicSpline<N-2>().Derivative(r) );
		}
	};

	template<>
	struct PolyharmonicSpline<1> {
		template<typename T>
		T operator()( T r ) const {
			return r;
		}
		template<typename T>
		T Derivative( T r ) const {
			return T(1);
		}
	};

	template<>
	struct PolyharmonicSpline<2> {
		template<typename T>
		T operator()( T r ) const {
			return r * r * log(r);
		}
		template<typename T>
		T Derivative( T r ) const {
			return r * ( T(1) + T(2) * log(r) );
		}
	};

	typedef PolyharmonicSpline<2> ThinPlateSpline;

	template<>
	struct PolyharmonicSpline<3> {
		template<typename T>
		T operator()( T r ) const {
			return r * r * r;
		}
		template<typename T>
		T Derivative( T r ) const {
			return T(3) * r * r;
		}
	};

	template<>
	struct PolyharmonicSpline<4> {
		template<typename T>
		T operator()( T r ) const {
			T r2 = r * r;
			return r2 * r2 * log(r);
		}
		template<typename T>
		T Derivative( T r ) const {
			T r2 = r * r;
			return r2 * r * ( T(1) + T(4) * log(r) );
		}
	};

	template<>
	struct PolyharmonicSpline<5> {
		template<typename T>
		T operator()( T r ) const {
			T r2 = r * r;
			return r2 * r2 * r;
		}
		template<typename T>
		T Derivative( T r ) const {
			T r2 = r * r;
			return T(5) * r2 * r2;
		}
	};

	template<>
	struct PolyharmonicSpline<6> {
		template<typename T>
		T operator()( T r ) const {
			T r3 = r * r * r;
			return r3 * r3;
		}
		template<typename T>
		T Derivative( T r ) const {
			T r2 = r * r;
			return r2 * r2 * r * ( T(1) + T(6) * log(r) );
		}
	};

	template<>
	struct PolyharmonicSpline<7> {
		template<typename T>
		T operator()( T r ) const {
			T r3 = r * r * r;
			return r3 * r3 * r;
		}
		template<typename T>
		T Derivative( T r ) const {
			T r3 = r * r * r;
			return T(7) * r3 * r3;
		}
	};

	template<>
	struct PolyharmonicSpline<8> {
		template<typename T>
		T operator()( T r ) const {
			T r2 = r * r;
			T r4 = r2 * r2;
			return r4 * r4;
		}
		template<typename T>
		T Derivative( T r ) const {
			T r3 = r * r * r;
			return r3 * r3 * r * ( T(1) + T(8) * log(r) );
		}
	};
}

#endif
