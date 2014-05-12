#ifndef INCLUDED_GEOMETRY_TRIG
#define INCLUDED_GEOMETRY_TRIG
#pragma warning (push)
#pragma warning (disable: 4985)
#include <cmath>
#pragma warning (pop)
#include "angle.h"

namespace DRDSP {

	template<typename T>
	inline T sin( Radians<T> x ) { return ::sin(x); }

	template<typename T>
	inline T cos( Radians<T> x ) { return ::cos(x); }

	template<typename T>
	inline T tan( Radians<T> x ) { return ::tan(x); }

	template<typename T>
	inline T csc( Radians<T> x ) { return T(1)/sin(x); }
	
	template<typename T>
	inline T sec( Radians<T> x ) { return T(1)/cos(x); }
	
	template<typename T>
	inline T cot( Radians<T> x ) { return T(1)/tan(x); }

	template<typename T>
	inline T csch( T x ) { return T(1)/sinh(x); }
	
	template<typename T>
	inline T sech( T x ) { return T(1)/cosh(x); }
	
	template<typename T>
	inline T coth( T x ) { return T(1)/tanh(x); }


	// derivatives

	template<typename T>
	inline T dsin( Radians<T> x ) { return cos(x); }
	
	template<typename T>
	inline T dcos( Radians<T> x ) { return -sin(x); }
	
	template<typename T>
	inline T dtan( Radians<T> x ) { T y = sec(x); return y*y; }

	template<typename T>
	inline T dcsc( Radians<T> x ) { return -csc(x)*cot(x); }
	
	template<typename T>
	inline T dsec( Radians<T> x ) { return sec(x)*tan(x); }
	
	template<typename T>
	inline T dcot( Radians<T> x ) { T y=csc(x); return -y*y; }

	template<typename T>
	inline T dasin( T x ) { return T(1)/sqrt(T(1)-x*x); }
	
	template<typename T>
	inline T dacos( T x ) { return -T(1)/sqrt(T(1)-x*x); }
	
	template<typename T>
	inline T datan( T x ) { return T(1)/(T(1)+x*x); }
	
	template<typename T>
	inline T dacsc( T x ) { return -T(1)/(abs(x)*sqrt(x*x-T(1))); }
	
	template<typename T>
	inline T dasec( T x ) { return T(1)/(abs(x)*sqrt(x*x-T(1))); }
	
	template<typename T>
	inline T dacot( T x ) { return -T(1)/(T(1)+x*x); }
	

	// hyperbolic derivatives

	template<typename T>
	inline T dsinh( T x ) { return cosh(x); }
	
	template<typename T>
	inline T dcosh( T x ) { return sinh(x); }
	
	template<typename T>
	inline T dtanh( T x ) { T y = sech(x); return y*y; }

	template<typename T>
	inline T dcsch( T x ) { return -csch(x)*coth(x); }
	
	template<typename T>
	inline T dsech( T x ) { return -sech(x)*tanh(x); }
	
	template<typename T>
	inline T dcoth( T x ) { T y=csch(x); return -y*y; }

	template<typename T>
	inline T dasinh( T x ) { return T(1)/sqrt(x*x+T(1)); }
	
	template<typename T>
	inline T dacosh( T x ) { return T(1)/sqrt(x*x-T(1)); }
	
	template<typename T>
	inline T datanh( T x ) { return T(1)/(T(1)-x*x); }
	
	template<typename T>
	inline T dacsch( T x ) { return -T(1)/(abs(x)*sqrt(T(1)+x*x)); }
	
	template<typename T>
	inline T dasech( T x ) { return -T(1)/(x*sqrt(T(1)-x*x)); }
	
	template<typename T>
	inline T dacoth( T x ) { return T(1)/(T(1)-x*x); }
	
}

#endif

