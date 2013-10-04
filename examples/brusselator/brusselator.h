#ifndef INCLUDED_BRUSSELATOR
#define INCLUDED_BRUSSELATOR
#include <DRDSP/types.h>
#include <DRDSP/dynamics/model_orig.h>

using namespace DRDSP;

struct Brusselator : ModelOriginal {
	double A, D1, D2;
	static const uint32_t nX = 32;
	static const uint32_t nY = 32;
	static const uint32_t N = nX * nY;
	double dx, dy;

	Brusselator();

	VectorXd VectorField( const VectorXd &X, const VectorXd &beta );
	MatrixXd VectorFieldD( const VectorXd &X, const VectorXd &beta );

	double VectorFieldCW( const VectorXd &x, const VectorXd &beta, uint32_t k );
	double VectorFieldDCW( const VectorXd &x, const VectorXd &beta, uint32_t k1, uint32_t k2 );

protected:
	inline double delta( uint32_t i, uint32_t j ) const;
	inline double X( const VectorXd &x, uint32_t i, uint32_t j ) const;
	inline double Y( const VectorXd &x, uint32_t i, uint32_t j ) const;
	double XDot( const VectorXd &x, const VectorXd &beta, uint32_t i, uint32_t j ) const;
	double YDot( const VectorXd &x, const VectorXd &beta, uint32_t i, uint32_t j ) const;
	double dXXDot( const VectorXd &x, const VectorXd &beta, uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const;
	double dYXDot( const VectorXd &x, uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const;
	double dXYDot( const VectorXd &x, const VectorXd &beta, uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const;
	double dYYDot( const VectorXd &x, uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const;
	double laplacianX( const VectorXd &x, uint32_t i, uint32_t j ) const;
	double laplacianY( const VectorXd &x, uint32_t i, uint32_t j ) const;
	double dAlaplacianA( uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const;

	uint32_t ip1( uint32_t i ) const;
	uint32_t im1( uint32_t i ) const;
	uint32_t jp1( uint32_t j ) const;
	uint32_t jm1( uint32_t j ) const;
};

#endif
