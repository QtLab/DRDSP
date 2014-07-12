#ifndef INCLUDED_FITZHUGH_NAGUMO
#define INCLUDED_FITZHUGH_NAGUMO
#include <DRDSP/dynamics/model.h>
#include <DRDSP/auto_diff.h>

using namespace DRDSP;

struct FitzHughNagumo : Model<> {
	int nX, nY, N;
	double sigma, lambda, kappa, tau, Du, Dv, dx, dy;
	complex<double> p, q;

	FitzHughNagumo() : FitzHughNagumo(32,32) {}

	FitzHughNagumo( int nX, int nY ) :
		Model<>(2*nX*nY),
		N(nX*nY), nX(nX), nY(nY), dx(1), dy(1),
		sigma(3), lambda(1), kappa(0), tau(3), Du(1), Dv(1) {}

	template<typename Derived>
	Matrix<typename Derived::Scalar,-1,1> operator()( const MatrixBase<Derived>& x ) const {
		typedef Derived::Scalar Scalar;
		auto u = compute_u(x);
		auto v = compute_v(x);
		Array<Scalar,-1,-1> uDot = action_potential(u) - sigma * v + Du * laplacian(u);
		Array<Scalar,-1,-1> vDot = ( u - v + Dv * laplacian(v) ) / tau;
		return vectorize( uDot, vDot );
	}

	SparseMatrix<double> Partials( const VectorXd& x ) const {
		return AutoDerivativeSparse( *this, x );
	}

protected:

	template<typename Derived>
	Array<typename Derived::Scalar,-1,-1> action_potential( const ArrayBase<Derived>& x ) const {
		return lambda * x - x * x * x - kappa;
	}

	template<typename Derived>
	Matrix<typename Derived::Scalar,-1,1> vectorize( const ArrayBase<Derived>& u, const ArrayBase<Derived>& v ) const {
		Matrix<Derived::Scalar,-1,1> x( stateDim );
		for(int j=0;j<nY;++j)
			x.segment(j*nX,nX) = u.col(j);
		for(int j=0;j<nY;++j)
			x.segment(N+j*nX,nX) = v.col(j);
		return x;
	}

	template<typename Derived>
	Array<typename Derived::Scalar,-1,-1> compute_u( const MatrixBase<Derived>& x ) const {
		Array<Derived::Scalar,-1,-1> u(nX,nY);
		for(int j=0;j<nY;++j)
			u.col(j) = x.segment(j*nX,nX);
		return u;
	}

	template<typename Derived>
	Array<typename Derived::Scalar,-1,-1> compute_v( const MatrixBase<Derived>& x ) const {
		Array<Derived::Scalar,-1,-1> v(nX,nY);
		for(int j=0;j<nY;++j)
			v.col(j) = x.segment(N+j*nX,nX);
		return v;
	}

	template<typename Derived>
	Array<typename Derived::Scalar,-1,-1> laplacian( const ArrayBase<Derived>& x ) const {
		typedef Derived::Scalar Scalar;
		Array<Scalar,-1,-1> LX(nX,nY), LY(nX,nY);
		
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

struct FitzHughNagumoFamily : Family<FitzHughNagumo> {
	int nX, nY;
	FitzHughNagumoFamily() : FitzHughNagumoFamily(32,32) {}
	FitzHughNagumoFamily( int nX, int nY ) : Family<FitzHughNagumo>(2*nX*nY,1), nX(nX), nY(nY) {}
	FitzHughNagumo operator()( const VectorXd& parameter ) const {
		FitzHughNagumo model(nX,nY);
		model.lambda = parameter[0];
		return model;
	}
};

#endif
