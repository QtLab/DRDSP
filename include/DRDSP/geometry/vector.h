#ifndef INCLUDED_GEOMETRY_VECTOR
#define INCLUDED_GEOMETRY_VECTOR
#include "../types.h"

namespace DRDSP {
	template<typename T> struct Traits;

	template<typename Derived>
	struct Vector {
		typedef typename Traits<Derived>::TScalar TScalar;
		virtual Derived operator+( const Derived &V ) const = 0;
		virtual Derived& operator+=( const Derived &V ) = 0;
		virtual Derived operator-( const Derived &V ) const = 0; 
		virtual Derived& operator-=( const Derived &V ) = 0;
		virtual Derived operator*( const TScalar &a ) const = 0;
		virtual Derived& operator*=( const TScalar &a ) = 0;
		virtual Derived operator/( const TScalar &a ) const = 0;
		virtual Derived& operator/=( const TScalar &a ) = 0;
		virtual Derived operator-() const = 0;
	};

	template<typename Derived>
	struct CoVector : Vector<Derived> {
		typedef typename Traits<Derived>::TVec TVec;
		virtual TScalar operator()( const TVec& V ) const = 0;
	};
	
	template<>
	struct Traits<double> {
		typedef double TScalar;
	};

	template<>
	struct Traits<float> {
		typedef float TScalar;
	};
	
	template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
	struct Traits<Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>> {
		typedef typename Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>::Scalar TScalar;
	};

	template<>
	struct Traits<struct CoVectorXd> {
		typedef VectorXd TVec;
		typedef Traits<TVec>::TScalar TScalar;
	};

	struct CoVectorXd : CoVector<CoVectorXd>, VectorXd {
		CoVectorXd( const VectorXd &rhs ) : CoVector<CoVectorXd>(), VectorXd(rhs) {}
		TScalar operator()( const TVec &V ) const {
			return dot(V);
		}
		CoVectorXd operator+( const CoVectorXd &M ) const {
			return VectorXd::operator+(M);
		}
		CoVectorXd& operator+=( const CoVectorXd &M ) {
			VectorXd::operator+=(M);
			return *this;
		}
		CoVectorXd operator-( const CoVectorXd &M ) const {
			return VectorXd::operator-(M);
		}
		CoVectorXd& operator-=( const CoVectorXd &M ) {
			VectorXd::operator-=(M);
			return *this;
		}
		CoVectorXd operator*( const TScalar &a ) const {
			return VectorXd::operator*(a);
		}
		CoVectorXd& operator*=( const TScalar &a ) {
			VectorXd::operator*=(a);
			return *this;
		}
		CoVectorXd operator/( const TScalar &a ) const {
			return VectorXd::operator/(a);
		}
		CoVectorXd& operator/=( const TScalar &a ) {
			VectorXd::operator/=(a);
			return *this;
		}
		CoVectorXd operator-() const {
			return VectorXd::operator-();
		}
	};

	template<>
	struct Traits<struct CoMatrixXd> {
		typedef MatrixXd TVec;
		typedef Traits<TVec>::TScalar TScalar;
	};

	struct CoMatrixXd : CoVector<CoMatrixXd>, MatrixXd {
		CoMatrixXd( const MatrixXd &rhs ) : CoVector<CoMatrixXd>(), MatrixXd(rhs) {}
		TScalar operator()( const TVec& V ) const {
			return (MatrixXd::operator*(V)).trace();
		}
		CoMatrixXd operator+( const CoMatrixXd &M ) const {
			return MatrixXd::operator+(M);
		}
		CoMatrixXd& operator+=( const CoMatrixXd &M ) {
			MatrixXd::operator+=(M);
			return *this;
		}
		CoMatrixXd operator-( const CoMatrixXd &M ) const {
			return MatrixXd::operator-(M);
		}
		CoMatrixXd& operator-=( const CoMatrixXd &M ) {
			MatrixXd::operator-=(M);
			return *this;
		}
		CoMatrixXd operator*( const TScalar &a ) const {
			return MatrixXd::operator*(a);
		}
		CoMatrixXd& operator*=( const TScalar &a ) {
			MatrixXd::operator*=(a);
			return *this;
		}
		CoMatrixXd operator/( const TScalar &a ) const {
			return MatrixXd::operator/(a);
		}
		CoMatrixXd& operator/=( const TScalar &a ) {
			MatrixXd::operator/=(a);
			return *this;
		}
		CoMatrixXd operator-() const {
			return MatrixXd::operator-();
		}
	};
}

#endif

