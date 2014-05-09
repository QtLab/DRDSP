#ifndef INCLUDED_EIGEN_DUAL
#define INCLUDED_EIGEN_DUAL
#include <limits>
#include <Eigen/Core>
#include "dual.h"

using namespace std;
using namespace Eigen;
using namespace DRDSP;

namespace Eigen {

	template<typename T>
	struct NumTraits<dual<T>> {
		enum {
			IsInteger = 0,
			IsSigned = 1,
			IsComplex = 0,
			RequireInitialization = NumTraits<T>::RequireInitialization,
			ReadCost = 2 * NumTraits<T>::ReadCost,
			AddCost = 2 * NumTraits<T>::AddCost,
			MulCost = 3 * NumTraits<T>::MulCost + NumTraits<T>::AddCost
		};
		typedef T Real;
		typedef dual<T> NonInteger;
		typedef dual<T> Nested;

		static inline Real epsilon() {
			return NumTraits<Real>::epsilon();
		}
		static inline Real dummy_precision() {
			return NumTraits<Real>::dummy_precision();
		}
		static inline T highest() {
			return dual<T>( NumTraits<T>::highest(), NumTraits<T>::highest() );
		}
		static inline T lowest() {
			return dual<T>( NumTraits<T>::lowest(), NumTraits<T>::lowest() );
		}
	};

	namespace internal {
		template<typename T> struct scalar_product_traits<T,dual<T>> {
			enum {
				// Cost = 2*NumTraits<T>::MulCost,
				Defined = 1
			};
			typedef dual<T> ReturnType;
		};

		template<typename T> struct scalar_product_traits<dual<T>,T> {
			enum {
				// Cost = 2*NumTraits<T>::MulCost,
				Defined = 1
			};
			typedef dual<T> ReturnType;
		};

		template<typename RealScalar> struct conj_helper<dual<RealScalar>, dual<RealScalar>, false,true> {
			typedef dual<RealScalar> Scalar;
			EIGEN_STRONG_INLINE Scalar pmadd(const Scalar& x, const Scalar& y, const Scalar& c) const
			{ return c + pmul(x,y); }
			EIGEN_STRONG_INLINE Scalar pmul(const Scalar& x, const Scalar& y) const
			{ return Scalar(numext::real(x)*numext::real(y) + numext::imag(x)*numext::imag(y), numext::imag(x)*numext::real(y) - numext::real(x)*numext::imag(y)); }
		};

		template<typename RealScalar> struct conj_helper<dual<RealScalar>, dual<RealScalar>, true,false> {
			typedef dual<RealScalar> Scalar;
			EIGEN_STRONG_INLINE Scalar pmadd(const Scalar& x, const Scalar& y, const Scalar& c) const
			{ return c + pmul(x,y); }
			EIGEN_STRONG_INLINE Scalar pmul(const Scalar& x, const Scalar& y) const
			{ return Scalar(numext::real(x)*numext::real(y) + numext::imag(x)*numext::imag(y), numext::real(x)*numext::imag(y) - numext::imag(x)*numext::real(y)); }
		};

		template<typename RealScalar> struct conj_helper<dual<RealScalar>, dual<RealScalar>, true,true> {
			typedef dual<RealScalar> Scalar;
			EIGEN_STRONG_INLINE Scalar pmadd(const Scalar& x, const Scalar& y, const Scalar& c) const
			{ return c + pmul(x,y); }
			EIGEN_STRONG_INLINE Scalar pmul(const Scalar& x, const Scalar& y) const
			{ return Scalar(numext::real(x)*numext::real(y) - numext::imag(x)*numext::imag(y), - numext::real(x)*numext::imag(y) - numext::imag(x)*numext::real(y)); }
		};

		template<typename RealScalar,bool Conj> struct conj_helper<dual<RealScalar>, RealScalar, Conj,false> {
			typedef dual<RealScalar> Scalar;
			EIGEN_STRONG_INLINE Scalar pmadd(const Scalar& x, const RealScalar& y, const Scalar& c) const
			{ return padd(c, pmul(x,y)); }
			EIGEN_STRONG_INLINE Scalar pmul(const Scalar& x, const RealScalar& y) const
			{ return conj_if<Conj>()(x)*y; }
		};

		template<typename RealScalar, bool Conj> struct conj_helper<RealScalar, dual<RealScalar>, false, Conj > {
			typedef dual<RealScalar> Scalar;
			EIGEN_STRONG_INLINE Scalar pmadd(const RealScalar& x, const Scalar& y, const Scalar& c) const
			{ return padd(c, pmul(x,y)); }
			EIGEN_STRONG_INLINE Scalar pmul(const RealScalar& x, const Scalar& y) const
			{ return x*conj_if<Conj>()(y); }
		};
	}
}

namespace DRDSP {

	template<typename T,int Rows,int Cols>
	using DualMatrix = Matrix<dual<T>,Rows,Cols>;

	typedef DualMatrix<float,-1,-1>  DualMatrixXf;
	typedef DualMatrix<double,-1,-1> DualMatrixXd;
	typedef DualMatrix<float,-1,1>  DualVectorXf;
	typedef DualMatrix<double,-1,1> DualVectorXd;

	template<typename Derived>
	DualMatrixXd Dualify( const MatrixBase<Derived>& x, const MatrixBase<Derived>& y ) {
		DualMatrixXd r( x.rows(), x.cols() );
		for(int64_t i=0;i<r.rows();++i) {
			for(int64_t j=0;j<r.cols();++j) {
				r(i,j).x = x(i,j);
				r(i,j).y = y(i,j);
			}
		}
		return r;
	}

	template<typename Derived>
	MatrixXd RealPart( const MatrixBase<Derived>& x ) {
		MatrixXd r( x.rows(), x.cols() );
		for(int64_t i=0;i<r.rows();++i) {
			for(int64_t j=0;j<r.cols();++j) {
				r(i,j) = x(i,j).x;
			}
		}
		return r;
	}

	template<typename Derived>
	MatrixXd DualPart( const MatrixBase<Derived>& x ) {
		MatrixXd r( x.rows(), x.cols() );
		for(int64_t i=0;i<r.rows();++i) {
			for(int64_t j=0;j<r.cols();++j) {
				r(i,j) = x(i,j).y;
			}
		}
		return r;
	}

}

#endif
