#ifndef INCLUDED_GINZBURG_LANDAU
#define INCLUDED_GINZBURG_LANDAU
#include <DRDSP/dynamics/model.h>
#include <DRDSP/auto_diff.h>

using namespace DRDSP;

struct GinzburgLandau : Model<> {
	int nX, nY, N;
	double gamma, dx, dy;
	complex<double> p, q;

	GinzburgLandau() : GinzburgLandau(32,32) {}

	GinzburgLandau( int nX, int nY ) : Model<>(2*nX*nY), N(nX*nY), nX(nX), nY(nY), p(0,1), q(0,1), gamma(1), dx(1), dy(1) {}

	template<typename Derived>
	Matrix<typename Derived::Scalar,-1,1> operator()( const MatrixBase<Derived>& x ) const {
		typedef Derived::Scalar Scalar;
		auto phi = order_parameter(x);
		return vectorize( Scalar(gamma) * phi + q * phi.abs2() * phi + complex<Scalar>(p) * laplacian(phi) );
	}

	SparseMatrix<double> Partials( const VectorXd& x ) const {
		return AutoDerivativeSparse( *this, x );
	}

protected:

	template<typename Derived>
	Matrix<typename Derived::Scalar::value_type,-1,1> vectorize( const ArrayBase<Derived>& phi ) const {
		Matrix<Derived::Scalar::value_type,-1,1> x( stateDim );
		for(int j=0;j<nY;++j)
			for(int i=0;i<nX;++i) {
				x[j*nX + i] = real( phi(i,j) );
				x[j*nX + i + 1] = imag( phi(i,j) );
			}
		return x;
	}

	template<typename Derived>
	Array<complex<typename Derived::Scalar>,-1,-1> order_parameter( const MatrixBase<Derived>& x ) const {
		Array<complex<Derived::Scalar>,-1,-1> phi(nX,nY);
		for(int j=0;j<nY;++j)
			for(int i=0;i<nX;++i) {
				phi(i,j).real( x[j*nX + i] );
				phi(i,j).imag( x[j*nX + i + 1] );
			}
		return phi;
	}

	template<typename Derived>
	Array<typename Derived::Scalar,-1,-1> laplacian( const ArrayBase<Derived>& x ) const {
		typedef Derived::Scalar Scalar;
		Array<Scalar,-1,-1> LX(nX,nY), LY(nX,nY);
		
		LX.row(0) = x.row(1) - Scalar(2.0) * x.row(0) + x.row(nX-1);
		for(int i=1;i<nX-1;++i)
			LX.row(i) = x.row(i+1) - Scalar(2.0) * x.row(i) + x.row(i-1);
		LX.row(nX-1) = x.row(0) - Scalar(2.0) * x.row(nX-1) + x.row(nX-2);

		LY.col(0) = x.col(1) - Scalar(2.0) * x.col(0) + x.col(nY-1);
		for(int j=1;j<nY-1;++j)
			LY.col(j) = x.col(j+1) - Scalar(2.0) * x.col(j) + x.col(j-1);
		LY.col(nY-1) = x.col(0) - Scalar(2.0) * x.col(nY-1) + x.col(nY-2);

		return LX * Scalar(1.0/(dx*dx)) + LY * Scalar(1.0/(dy*dy));
	}

};

struct GinzburgLandauFamily : Family<GinzburgLandau> {
	int nX, nY;
	GinzburgLandauFamily() : GinzburgLandauFamily(32,32) {}
	GinzburgLandauFamily( int nX, int nY ) : Family<GinzburgLandau>(2*nX*nY,1), nX(nX), nY(nY) {}
	GinzburgLandau operator()( const VectorXd& parameter ) const {
		GinzburgLandau model(nX,nY);
		model.gamma = parameter[0];
		return model;
	}
};

#endif
