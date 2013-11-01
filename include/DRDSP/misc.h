#ifndef INCLUDED_MISC
#define INCLUDED_MISC
#include "types.h"

namespace DRDSP {

	template<typename T>
	T Clamp( T x, T lower, T upper ) {
		if( x > upper ) return upper;
		if( x < lower ) return lower;
		return x;
	}

	template<typename T>
	void Wrap( T& x, T a, T b ) {
		T length = b - a;
		while( x < a ) x += length;
		while( x > b ) x -= length;
	}

	template<typename T>
	bool IsEven( T x ) {
		return (x % 2) == 0;
	}

	template<typename T>
	bool IsOdd( T x ) {
		return (x % 2) != 0;
	}

	template<typename T>
	T Delta( uint32_t i, uint32_t j ) {
		return (i==j)?T(1):T(0);
	}

}

#endif
