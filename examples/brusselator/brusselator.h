#ifndef INCLUDED_BRUSSELATOR
#define INCLUDED_BRUSSELATOR
#include <DRDSP/types.h>
#include <DRDSP/dynamics/model.h>

using namespace DRDSP;

struct Brusselator : Model<> {
	double A, B, D1, D2;
	static const uint32_t nX = 32;
	static const uint32_t nY = 32;
	static const uint32_t N = nX * nY;
	double dx, dy;

	Brusselator() : Brusselator(1.0) {}

	explicit Brusselator( double B ) : Model<>(2*N), A(1), B(B), D1(1), D2(1), dx(1), dy(1) {}

	VectorXd operator()( const VectorXd &X ) const;
	MatrixXd Partials( const VectorXd &X ) const;

	double operator()( const VectorXd &x, uint32_t k ) const;
	double Partials( const VectorXd &x, uint32_t k1, uint32_t k2 ) const;

protected:
	inline double delta( uint32_t i, uint32_t j ) const;
	inline double X( const VectorXd& x, uint32_t i, uint32_t j ) const;
	inline double Y( const VectorXd& x, uint32_t i, uint32_t j ) const;
	double XDot( const VectorXd& x, uint32_t i, uint32_t j ) const;
	double YDot( const VectorXd& x, uint32_t i, uint32_t j ) const;
	double dXXDot( const VectorXd& x, uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const;
	double dYXDot( const VectorXd& x, uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const;
	double dXYDot( const VectorXd& x, uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const;
	double dYYDot( const VectorXd& x, uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const;
	double laplacianX( const VectorXd& x, uint32_t i, uint32_t j ) const;
	double laplacianY( const VectorXd& x, uint32_t i, uint32_t j ) const;
	double dAlaplacianA( uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const;

	uint32_t ip1( uint32_t i ) const;
	uint32_t im1( uint32_t i ) const;
	uint32_t jp1( uint32_t j ) const;
	uint32_t jm1( uint32_t j ) const;
};

struct BrusselatorFamily : Family<Brusselator> {
	BrusselatorFamily() : Family<Brusselator>(2*Brusselator::N,1) {}
	Brusselator operator()( const VectorXd& parameter ) const {
		return Brusselator( parameter[0] );
	}
};

#endif
