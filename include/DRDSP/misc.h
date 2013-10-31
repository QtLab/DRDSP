#ifndef INCLUDED_MISC
#define INCLUDED_MISC

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

}

#endif
