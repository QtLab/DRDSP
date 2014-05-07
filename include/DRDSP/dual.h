#ifndef INCLUDED_DUAL
#define INCLUDED_DUAL
#include <cmath>

namespace DRDSP {

	template<typename T>
	struct dual {
		T x, y;

		dual() = default;

		explicit dual( T real ) : x(real), y(0) {}

		dual( T real, T dual ) : x(real), y(dual) {}

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

		dual<T> operator*( T a ) const {
			return dual<T>(x*a,y*a);
		}

		dual<T>& operator*=( T a ) {
			x *= a;
			y *= a;
			return *this;
		}

		dual<T> operator/( T a ) const {
			return dual<T>(x/a,y/a);
		}

		dual<T>& operator/=( T a ) {
			x /= a;
			y /= a;
			return *this;
		}

		T real_part() const {
			return x;
		}

		T dual_part() const {
			return y;
		}
	};

	template<typename T>
	dual<T> operator*( T a, const dual<T>& rhs ) {
		return rhs*a;
	}

	template<typename T>
	dual<T> conjugate( const dual<T>& d ) {
		return dual<T>(d.x,-d.y);
	}

	template<typename T>
	T modulus( const dual<T>& d ) {
		return d.x;
	}

	template<typename T>
	T arg( const dual<T>& d ) {
		return d.y/d.x;
	}

	template<typename T>
	dual<T> sqrt( const dual<T>& d ) {
		T a = sqrt(d.x);
		return dual<T>(a,d.y/(a*T(2)));
	}

	template<typename T>
	dual<T> exp( const dual<T>& d ) {
		T ex = exp(d.x);
		return dual<T>(ex,ex*d.y);
	}

	template<typename T>
	dual<T> log( const dual<T>& d ) {
		T m = modulus(d);
		T a = Arg(d);
		return dual<T>(log(m),a);
	}

	template<typename T>
	dual<T> pow( const dual<T>& d, const dual<T>& p ) {
		T r = modulus(d);
		T rp = pow(r,p.x);
		return dual<T>(rp,rp * (p.x*Arg(d)+p.y*log(r)) );
	}

	template<typename T>
	dual<T> pow( const dual<T>& d, T p ) {
		if( d.x == T(0) ) {
			return dual<T>(0,0);
		}
		T rp = pow(d.x,p);
		return dual<T>(rp,rp*p*d.y/d.x);
	}

	template<typename T>
	dual<T> sin( const dual<T> &d ) {
		return dual<T>(sin(d.x),d.y*cos(d.x));
	}

	template<typename T>
	dual<T> cos( const dual<T>& d ) {
		return dual<T>(cos(d.x),-d.y*sin(d.x));
	}

	template<typename T>
	dual<T> Tan( const dual<T>& d ) {
		T t = tan(d.x);
		return dual<T>(t,d.y*(T(1)+t*t));
	}

	template<typename T>
	dual<T> sec( const dual<T>& d ) {
		T s = sec(d.x);
		return dual<T>(s,d.y*s*tan(d.x));
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
		return dual<T>( asin(d.x), d.y*dasin(x) );
	}

	template<typename T>
	dual<T> acos( const dual<T>& d ) {
		return dual<T>( acos(d.x), d.y*dacos(x) );
	}

	template<typename T>
	dual<T> atan( const dual<T>& d ) {
		return dual<T>( atan(d.x), d.y*datan(x) );
	}

	template<typename T>
	dual<T> asec( const dual<T>& d ) {
		return dual<T>( asec(d.x), d.y*dasec(x) );
	}

	template<typename T>
	dual<T> acsc( const dual<T>& d ) {
		return dual<T>( acsc(d.x), d.y*dacsc(x) );
	}

	template<typename T>
	dual<T> acot( const dual<T>& d ) {
		return dual<T>( acot(d.x), d.y*dacot(x) );
	}

	typedef dual<double> duald;
	typedef dual<float>  dualf;

}

#endif
