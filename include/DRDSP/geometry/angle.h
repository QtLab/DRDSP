#ifndef INCLUDED_GEOMETRY_ANGLE
#define INCLUDED_GEOMETRY_ANGLE

namespace DRDSP {

	const double PI = 3.1415926535897932384;
	const double TAU = 2.0 * PI;

	const double RAD_PER_DEG = (PI/180.0);
	const double DEG_PER_RAD = (180.0/PI);
	
	template<typename T>
	struct Radians;

	template<typename T>
	struct Degrees {
		T value;
		
		Degrees() = default;
		
		Degrees( T rhs ) : value(rhs) {}
		
		Degrees( const Radians<T>& rhs ) : value( rhs.value * T(DEG_PER_RAD) ) {}

		Degrees<T>& operator=( const Radians<T>& rhs ) {
			value = rhs.value * T(DEG_PER_RAD);
			return *this;
		}
		
		Degrees<T> operator+( const Degrees<T>& rhs ) const {
			return Degrees<T>( value + rhs.value );
		}
		
		Degrees<T> operator-( const Degrees<T>& rhs ) const {
			return Degrees<T>( value - rhs.value );
		}
		
		Degrees<T> operator*( T x ) const {
			return Degrees<T>( value * x );
		}
		
		Degrees<T> operator/( T x ) const {
			return Degrees<T>( value / x );
		}
		
		Degrees<T>& operator+=( const Degrees<T>& rhs ) {
			value += rhs.value;
			return *this;
		}
		
		Degrees<T>& operator-=( const Degrees<T>& rhs ) {
			value -= rhs.value;
			return *this;
		}
		
		Degrees<T>& operator*=( T x ) {
			value *= x;
			return *this;
		}
		
		Degrees<T>& operator/=( T x ) {
			value /= x;
			return *this;
		}
		
		Degrees<T> operator-() const {
			return Degrees<T>( -value );
		}

		bool operator>( const Degrees<T>& rhs ) const {
			return value > rhs.value;
		}

		bool operator<( const Degrees<T>& rhs ) const {
			return value < rhs.value;
		}

		bool operator>=( const Degrees<T>& rhs ) const {
			return value >= rhs.value;
		}

		bool operator<=( const Degrees<T>& rhs ) const {
			return value <= rhs.value;
		}

		operator T() const {
			return value;
		}

		/*Degrees<T> operator"" deg( T x ) const {
			return Degrees<T>( x );
		}*/

	};

	typedef Degrees<float>  Degreesf;
	typedef Degrees<double> Degreesd;

	template<typename T>
	struct Radians {
		T value;
		
		Radians() = default;
		
		Radians( T rhs ) : value(rhs) {}
		
		Radians( const Degrees<T>& rhs ) : value( rhs.value * T(RAD_PER_DEG) ) {}

		Radians<T>& operator=( const Degrees<T>& rhs ) {
			value = rhs.value * T(RAD_PER_DEG);
			return *this;
		}
		
		Radians<T> operator+( const Radians<T>& rhs ) const {
			return Radians<T>( value + rhs.value );
		}
		
		Radians<T> operator-( const Radians<T>& rhs ) const {
			return Radians<T>( value - rhs.value );
		}
		
		Radians<T> operator*( T x ) const {
			return Radians<T>( value * x );
		}
		
		Radians<T> operator/( T x ) const {
			return Radians<T>( value / x );
		}
		
		Radians<T>& operator+=( const Radians<T>& rhs ) {
			value += rhs.value;
			return *this;
		}
		
		Radians<T>& operator-=( const Radians<T>& rhs ) {
			value -= rhs.value;
			return *this;
		}
		
		Radians<T>& operator*=( T x ) {
			value *= x;
			return *this;
		}
		
		Radians<T>& operator/=( T x ) {
			value /= x;
			return *this;
		}
		
		Radians<T> operator-() const {
			return Radians<T>( -value );
		}

		bool operator>( const Radians<T>& rhs ) const {
			return value > rhs.value;
		}

		bool operator<( const Radians<T>& rhs ) const {
			return value < rhs.value;
		}

		bool operator>=( const Radians<T>& rhs ) const {
			return value >= rhs.value;
		}

		bool operator<=( const Radians<T>& rhs ) const {
			return value <= rhs.value;
		}

		operator T() const {
			return value;
		}

		/*Radians<T> operator"" rad( T x ) const {
			return Radians<T>( x );
		}*/

	};

	typedef Radians<float>  Radiansf;
	typedef Radians<double> Radiansd;

}

#endif
