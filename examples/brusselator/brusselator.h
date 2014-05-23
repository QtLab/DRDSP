#ifndef INCLUDED_BRUSSELATOR
#define INCLUDED_BRUSSELATOR
#include <DRDSP/types.h>
#include <DRDSP/dynamics/model.h>
#include <DRDSP/auto_diff.h>

using namespace DRDSP;

struct Brusselator : Model<> {
	static const int nX = 32;
	static const int nY = 32;
	static const int N = nX * nY;

	double A, B, D1, D2, dx, dy;

	Brusselator() : Brusselator(1.0) {}

	explicit Brusselator( double B ) : Model<>(2*N), A(1), B(B), D1(1), D2(1), dx(1), dy(1) {}

	template<typename Derived>
	Matrix<typename Derived::Scalar,-1,1> operator()( const MatrixBase<Derived>& x ) const {
		typedef Derived::Scalar Scalar;
		auto X = compute_X(x);
		auto Y = compute_Y(x);
		Array<Scalar,-1,-1> XXY = X*X*Y;
		Array<Scalar,-1,-1> XDot = A + XXY - (B+1.0)*X + D1 * laplacian(X);
		Array<Scalar,-1,-1> YDot = B*X - XXY + D2 * laplacian(Y);

		Matrix<Scalar,-1,1> r(dimension);

		for(int j=0;j<nY;++j)
			r.segment(j*nX,nX) = XDot.col(j);

		for(int j=0;j<nY;++j)
			r.segment(N+j*nX,nX) = YDot.col(j);

		return r;
	}

	SparseMatrix<double> Partials( const VectorXd& x ) const {
		return AutoDerivativeSparse( *this, x );
	}

protected:

	template<typename Derived>
	Array<typename Derived::Scalar,-1,-1> compute_X( const MatrixBase<Derived>& x ) const {
		Array<Derived::Scalar,-1,-1> X(nX,nY);
		for(int j=0;j<nY;++j)
			X.col(j) = x.segment(j*nX,nX);
		return X;
	}

	template<typename Derived>
	Array<typename Derived::Scalar,-1,-1> compute_Y( const MatrixBase<Derived>& x ) const {
		Array<Derived::Scalar,-1,-1> Y(nX,nY);
		for(int j=0;j<nY;++j)
			Y.col(j) = x.segment(N+j*nX,nX);
		return Y;
	}

	template<typename Derived>
	Array<typename Derived::Scalar,-1,-1> laplacian( const ArrayBase<Derived>& x ) const {
		Array<Derived::Scalar,-1,-1> LX(nX,nY), LY(nX,nY);
		
		LX.row(0) = x.row(1) - 2.0 * x.row(0) + x.row(nX-1);
		for(int i=1;i<nX-1;++i)
			LX.row(i) = x.row(i+1) - 2.0 * x.row(i) + x.row(i-1);
		LX.row(nX-1) = x.row(0) - 2.0 * x.row(nX-1) + x.row(nX-2);

		LY.col(0) = x.col(1) - 2.0 * x.col(0) + x.col(nY-1);
		for(int j=1;j<nY-1;++j)
			LY.col(j) = x.col(j+1) - 2.0 * x.col(j) + x.col(j-1);
		LY.col(nY-1) = x.col(0) - 2.0 * x.col(nY-1) + x.col(nY-2);

		return LX * (1.0/(dx*dx)) + LY * (1.0/(dy*dy));
	}

};

struct BrusselatorFamily : Family<Brusselator> {
	BrusselatorFamily() : Family<Brusselator>(2*Brusselator::N,1) {}
	Brusselator operator()( const VectorXd& parameter ) const {
		return Brusselator( parameter[0] );
	}
};

#endif
