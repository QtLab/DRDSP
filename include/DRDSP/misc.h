#ifndef INCLUDED_MISC
#define INCLUDED_MISC
#include "types.h"
#include <numeric>
#include <algorithm>
#include <vector>

namespace DRDSP {

	std::vector<VectorXd> ParameterList( const VectorXd& pMin, const VectorXd& pMax, uint32_t N );

	std::vector<VectorXd> ParameterList( double pMin, double pMax, uint32_t N );

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
	inline bool IsEven( T x ) {
		return (x % 2) == 0;
	}

	template<typename T>
	inline bool IsOdd( T x ) {
		return (x % 2) != 0;
	}

	template<typename T>
	inline T Delta( uint32_t i, uint32_t j ) {
		return (i==j)?T(1):T(0);
	}

	template<typename C,typename T>
	inline T accumulate( const C& container, T&& value ) {
		return std::accumulate(
			std::cbegin( container ),
			std::cend( container ),
			std::forward<T>(value)
		);
	}

	template<typename C>
	inline auto minmax_element( const C& container ) -> decltype(std::minmax_element(std::cbegin(container),std::cend(container))) {
		return std::minmax_element(
			std::cbegin( container ),
			std::cend( container )
		);
	}

	template<typename Derived>
	Matrix<typename Derived::Scalar,-1,1> Vectorize( const MatrixBase<Derived>& A ) {
		Matrix<typename Derived::Scalar,-1,1> V( A.size() );
		for(int64_t i=0;i<A.cols();++i) {
			V.segment(i*A.rows(),A.rows()) = A.col(i);
		}
		return V;
	}

}

#endif
