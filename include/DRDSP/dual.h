#ifndef INCLUDED_DUAL
#define INCLUDED_DUAL
/**
 * An implementation of dual numbers that has the same structure as std::complex.
 *
 */
#include "geometry/trig.h"
#include <complex>

namespace DRDSP {

	template<typename T>
	struct dual {
		T x, y;

		dual( T re = T(), T du = T() ) : x(re), y(du) {}

		dual<T> operator-() const {
			return dual<T>(-x,-y);
		}

		dual<T> operator+( const dual<T>& rhs ) const {
			return dual<T>(x+rhs.x,y+rhs.y);
		}

		dual<T> operator-( const dual<T>& rhs ) const {
			return dual<T>(x-rhs.x,y-rhs.y);
		}

		dual<T> operator*( const dual<T>& rhs ) const {
			dual<T> res;
			res.x = x*rhs.x;
			res.y = y*rhs.x + x*rhs.y;
			return res;
		}

		dual<T> operator/( const dual<T>& rhs ) const {
			dual<T> res;
			res.x = x/rhs.x;
			res.y = (y*rhs.x - x*rhs.y)/(rhs.x*rhs.x);
			return res;
		}

		dual<T>& operator+=( const dual<T>& rhs ) {
			x += rhs.x;
			y += rhs.y;
			return *this;
		}

		dual<T>& operator-=( const dual<T>& rhs ) {
			x -= rhs.x;
			y -= rhs.y;
			return *this;
		}

		dual<T>& operator*=( const dual<T>& rhs ) {
			*this = *this * rhs;
			return *this;
		}

		dual<T>& operator/=( const dual<T>& rhs ) {
			*this = *this / rhs;
			return *this;
		}

		dual<T>& operator=( T rhs ) {
			x = rhs;
			y = T(0);
			return *this;
		}

		dual<T> operator+( T rhs ) const {
			return dual<T>(x+rhs,y);
		}

		dual<T> operator-( T rhs ) const {
			return dual<T>(x-rhs,y);
		}

		dual<T> operator*( T rhs ) const {
			return dual<T>(x*rhs,y*rhs);
		}

		dual<T> operator/( T rhs ) const {
			return dual<T>(x/rhs,y/rhs);
		}
		
		dual<T>& operator+=( T rhs ) {
			x += rhs;
			return *this;
		}
		
		dual<T>& operator-=( T rhs ) {
			x -= rhs;
			return *this;
		}

		dual<T>& operator*=( T rhs ) {
			x *= rhs;
			y *= rhs;
			return *this;
		}

		dual<T>& operator/=( T rhs ) {
			x /= rhs;
			y /= rhs;
			return *this;
		}

		T real() const {
			return x;
		}

		T dual_part() const {
			return y;
		}

		void real( T re ) {
			x = re;
		}

		void dual_part( T du ) {
			y = du;
		}

		bool operator==( const dual<T>& rhs ) const {
			return x == rhs.x && y == rhs.y;
		}

		bool operator!=( const dual<T>& rhs ) const {
			return !(*this == rhs);
		}

	};

	template<typename T>
	dual<T> operator+( T lhs, const dual<T>& rhs ) {
		return rhs + lhs;
	}

	template<typename T>
	dual<T> operator-( T lhs, const dual<T>& rhs ) {
		return -rhs + lhs;
	}

	template<typename T>
	dual<T> operator*( T lhs, const dual<T>& rhs ) {
		return rhs * lhs;
	}

	template<typename T>
	dual<T> operator/( T lhs, const dual<T>& rhs ) {
		return dual<T>(lhs) / rhs;
	}

	template<typename T>
	T real( const dual<T>& d ) {
		return d.real();
	}

	template<typename T>
	T dual_part( const dual<T>& d ) {
		return d.dual_part();
	}

	template<typename T>
	dual<T> abs( const dual<T>& d ) {
		T m = abs(d.x);
		return dual<T>(m,d.y*d.x/m);
	}

	template<typename T>
	T arg( const dual<T>& d ) {
		return d.y/d.x;
	}

	template<typename T>
	dual<T> sqrt( const dual<T>& d ) {
		T a = ::sqrt(d.x);
		return dual<T>(a,d.y/(a*T(2)));
	}

	template<typename T>
	dual<T> exp( const dual<T>& d ) {
		T ex = ::exp(d.x);
		return dual<T>(ex,ex*d.y);
	}

	template<typename T>
	dual<T> log( const dual<T>& d ) {
		return dual<T>(::log(d.x),d.y/d.x);
	}

	template<typename T>
	dual<T> pow( const dual<T>& d, const dual<T>& p ) {
		if( d.x == T(0) ) return dual<T>(0,0);
		T rp = ::pow(d.x,p.x);
		return dual<T>( rp, rp*( p.y*::log(d.x) + p.x*d.y/d.x ) );
	}

	template<typename T>
	dual<T> pow( const dual<T>& d, T p ) {
		if( d.x == T(0) ) return dual<T>(0,0);
		T rp = ::pow(d.x,p.x);
		return dual<T>( rp, (rp*p*d.y)/d.x );
	}

	template<typename T>
	dual<T> sin( const dual<T> &d ) {
		return dual<T>(::sin(d.x),d.y*::cos(d.x));
	}

	template<typename T>
	dual<T> cos( const dual<T>& d ) {
		return dual<T>(::cos(d.x),-d.y*::sin(d.x));
	}

	template<typename T>
	dual<T> tan( const dual<T>& d ) {
		T t = ::tan(d.x);
		return dual<T>(t,d.y*(T(1)+t*t));
	}

	template<typename T>
	dual<T> sec( const dual<T>& d ) {
		T s = sec(d.x);
		return dual<T>(s,d.y*s*::tan(d.x));
	}

	template<typename T>
	dual<T> csc( const dual<T>& d ) {
		T s = csc(d.x);
		return dual<T>(s,-d.y*s*cot(d.x));
	}

	template<typename T>
	dual<T> cot( const dual<T>& d ) {
		T t = cot(d.x);
		return dual<T>(t,-d.y*(T(1)+t*t));
	}

	template<typename T>
	dual<T> asin( const dual<T>& d ) {
		return dual<T>( ::asin(d.x), d.y*dasin(d.x) );
	}

	template<typename T>
	dual<T> acos( const dual<T>& d ) {
		return dual<T>( ::acos(d.x), d.y*dacos(d.x) );
	}

	template<typename T>
	dual<T> atan( const dual<T>& d ) {
		return dual<T>( ::atan(d.x), d.y*datan(d.x) );
	}

	template<typename T>
	dual<T> asec( const dual<T>& d ) {
		return dual<T>( asec(d.x), d.y*dasec(d.x) );
	}

	template<typename T>
	dual<T> acsc( const dual<T>& d ) {
		return dual<T>( acsc(d.x), d.y*dacsc(d.x) );
	}

	template<typename T>
	dual<T> acot( const dual<T>& d ) {
		return dual<T>( acot(d.x), d.y*dacot(d.x) );
	}

	template<typename T>
	dual<T> sinh( const dual<T>& d ) {
		return dual<T>( ::sinh(d.x), d.y*::cosh(d.x) );
	}

	template<typename T>
	dual<T> cosh( const dual<T>& d ) {
		return dual<T>( ::cosh(d.x), d.y*::sinh(d.x) );
	}

	template<typename T>
	dual<T> tanh( const dual<T>& d ) {
		T t = ::tanh(d.x);
		return dual<T>( t, d.y*(T(1)-t*t) );
	}

	template<typename T>
	dual<T> coth( const dual<T>& d ) {
		T t = coth(d.x);
		return dual<T>( t, d.y*(T(1)-t*t) );
	}

	template<typename T>
	dual<T> csch( const dual<T>& d ) {
		T c = csch(d.x);
		return dual<T>( c, -d.y*c*coth(d.x) );
	}

	template<typename T>
	dual<T> sech( const dual<T>& d ) {
		T s = sech(d.x);
		return dual<T>( s, -d.y*s*::tanh(d.x) );
	}

	template<typename T>
	dual<T> asinh( const dual<T>& d ) {
		return dual<T>( ::asinh(d.x), d.y*dasinh(d.x) );
	}

	template<typename T>
	dual<T> acosh( const dual<T>& d ) {
		return dual<T>( ::acosh(d.x), d.y*dacosh(d.x) );
	}

	template<typename T>
	dual<T> atanh( const dual<T>& d ) {
		return dual<T>( ::atanh(d.x), d.y*datanh(d.x) );
	}

	template<typename T>
	dual<T> acsch( const dual<T>& d ) {
		return dual<T>( acsch(d.x), d.y*dacsch(d.x) );
	}

	template<typename T>
	dual<T> asech( const dual<T>& d ) {
		return dual<T>( asech(d.x), d.y*dasech(d.x) );
	}

	template<typename T>
	dual<T> acoth( const dual<T>& d ) {
		return dual<T>( acoth(d.x), d.y*dacoth(d.x) );
	}

	typedef dual<double> duald;
	typedef dual<float>  dualf;

}

#endif
