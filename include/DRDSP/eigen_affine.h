#ifndef INCLUDED_EIGEN_AFFINE
#define INCLUDED_EIGEN_AFFINE
#include <Eigen/Core>

using namespace Eigen;

namespace DRDSP {

	template<typename Scalar, int Rows, int Cols>
	struct Affine {
		typedef Scalar Scalar;
		enum {
			Rows = Rows,
			Cols = Cols
		};

		uint32_t sourceDim, targetDim;
		Matrix<Scalar,Rows,Cols> linear;
		Matrix<Scalar,Rows,1> translation;

		Affine() = default;

		Affine( int64_t rows, int64_t cols ) :
			linear(rows,cols),
			translation(rows),
			sourceDim((uint32_t)cols),
			targetDim((uint32_t)rows)
		{}

		template<typename Derived>
		Matrix<Scalar,Rows,1> operator()( const MatrixBase<Derived>& x ) const {
			return linear * x + translation;
		}

		Affine<Scalar,Rows,Cols>& setIdentity() {
			linear.setIdentity();
			translation.setZero();
			return *this;
		}

		Affine<Scalar,Rows,Cols>& setZero() {
			linear.setZero();
			translation.setZero();
			return *this;
		}

		Affine<Scalar,Rows,Cols>& setIdentity( int64_t rows, int64_t cols ) {
			linear.setIdentity(rows,cols);
			translation.setZero(rows);
			return *this;
		}

		Affine<Scalar,Rows,Cols>& setZero( int64_t rows, int64_t cols ) {
			linear.setZero(rows,cols);
			translation.setZero(rows);
			return *this;
		}

	};
	
	typedef Affine<float, -1,-1> AffineXf;
	typedef Affine<double,-1,-1> AffineXd;

	AffineXd VecToAffine( const VectorXd& param, int64_t rows ) {
		int64_t cols = param.size() / rows - 1;
		AffineXd A( rows, cols );
		for(int64_t i=0;i<cols;++i) {
			A.linear.col(i) = param.segment(i*rows,rows);
		}
		A.translation = param.segment(cols*rows,rows);
		return A;
	}

	VectorXd AffineToVec( const AffineXd& A ) {
		int64_t rows = A.linear.rows();
		int64_t cols = A.linear.cols();
		VectorXd v( rows * ( cols + 1 ) );
		for(int64_t i=0;i<cols;++i) {
			v.segment(i*rows,rows) = A.linear.col(i);
		}
		v.segment(cols*rows,rows) = A.translation;
		return v;
	}

}

#endif
