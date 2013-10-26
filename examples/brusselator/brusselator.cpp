#include "brusselator.h"

Brusselator::Brusselator() : ModelParameterized(2*N,1) {}

VectorXd Brusselator::VectorField( const VectorXd &x, const VectorXd &beta ) {
	VectorXd r(dimension);
	uint32_t k=0;
	for(uint32_t j=0;j<nY;j++)
		for(uint32_t i=0;i<nX;i++) {
			r(k++) = XDot(x,beta,i,j);
			r(k++) = YDot(x,beta,i,j);
		}
	return r;
}

MatrixXd Brusselator::Partials( const VectorXd &x, const VectorXd &beta ) {
	MatrixXd r((int)dimension,(int)dimension);
	uint32_t m=0,n=0;
	for(uint32_t j=0;j<nY;j++)
		for(uint32_t i=0;i<nX;i++) {
			n = 0;
			for(uint32_t nu=0;nu<nY;nu++)
				for(uint32_t mu=0;mu<nX;mu++) {
					r(m,n++) = dXXDot(x,beta,i,j,mu,nu);
					r(m,n++) = dYXDot(x,i,j,mu,nu);
				}
			m++;
			n = 0;
			for(uint32_t nu=0;nu<nY;nu++)
				for(uint32_t mu=0;mu<nX;mu++) {
					r(m,n++) = dXYDot(x,beta,i,j,mu,nu);
					r(m,n++) = dYYDot(x,i,j,mu,nu);
				}
			m++;
		}
	return r;
}

BrusselatorCW::BrusselatorCW() : ModelParameterizedCW(2*N,1) {}

double BrusselatorCW::VectorField( const VectorXd &x, const VectorXd &beta, uint32_t k ) {
	uint32_t i,j;
	if( !(k % 2) ) {
		j = k / (2*nX);
		i = k/2 - j * nX;
		return XDot(x,beta,i,j);
	} else {
		k -= 1;
		j = k / (2*nX);
		i = k/2 - j * nX;
		return YDot(x,beta,i,j);
	}
}

double BrusselatorCW::Partials( const VectorXd &x, const VectorXd &beta, uint32_t k1, uint32_t k2 ) {
	uint32_t i1,j1,i2,j2;
	if( !(k1 % 2) ) {
		j1 = k1 / (2*nX);
		i1 = k1/2 - j1 * nX;
		if( !(k2 % 2) ) {
			j2 = k2 / (2*nX);
			i2 = k2/2 - j2 *nX;
			return dXXDot(x,beta,i1,j1,i2,j2);
		} else {
			k2 -= 1;
			j2 = k2 / (2*nX);
			i2 = k2/2 - j2 * nX;
			return dYXDot(x,i1,j1,i2,j2);
		}
	} else {
		k1 -= 1;
		j1 = k1 / (2*nX);
		i1 = k1/2 - j1 * nX;
		if( !(k2 % 2) ) {
			j2 = k2 / (2*nX);
			i2 = k2/2 - j2 * nX;
			return dXYDot(x,beta,i1,j1,i2,j2);
		} else {
			k2 -= 1;
			j2 = k2 / (2*nX);
			i2 = k2/2 - j2 * nX;
			return dYYDot(x,i1,j1,i2,j2);
		}
	}
}

BrusselatorBase::BrusselatorBase() : A(1), D1(1), D2(1), dx(1), dy(1) {}

inline double BrusselatorBase::delta( uint32_t i, uint32_t j ) const {
	return (i==j)?1.0:0.0;
}

inline double BrusselatorBase::X( const VectorXd &x, uint32_t i, uint32_t j ) const {
	return x(2*(nX*j+i));
}

inline double BrusselatorBase::Y( const VectorXd &x, uint32_t i, uint32_t j ) const {
	return x(2*(nX*j+i)+1);
}

double BrusselatorBase::XDot( const VectorXd &x, const VectorXd &beta, uint32_t i, uint32_t j ) const {
	double B = beta(0);
	return A + X(x,i,j)*X(x,i,j)*Y(x,i,j) - (B+1.0)*X(x,i,j) + D1 * laplacianX(x,i,j);
}

double BrusselatorBase::YDot( const VectorXd &x, const VectorXd &beta, uint32_t i, uint32_t j ) const {
	double B = beta(0);
	return B*X(x,i,j) - X(x,i,j)*X(x,i,j)*Y(x,i,j) + D2 * laplacianY(x,i,j);
}

double BrusselatorBase::dXXDot( const VectorXd &x, const VectorXd &beta, uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const {
	double B = beta(0);
	return (2.0*X(x,i,j)*Y(x,i,j)-(B+1.0))*delta(i,mu)*delta(j,nu) + D1 * dAlaplacianA(i,j,mu,nu);
}

double BrusselatorBase::dYXDot( const VectorXd &x, uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const {
	return X(x,i,j)*X(x,i,j)*delta(i,mu)*delta(j,nu);
}

double BrusselatorBase::dXYDot( const VectorXd &x, const VectorXd &beta, uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const {
	double B = beta(0);
	return (B - 2.0*X(x,i,j)*Y(x,i,j))*delta(i,mu)*delta(j,nu);
}

double BrusselatorBase::dYYDot( const VectorXd &x, uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const {
	return -X(x,i,j)*X(x,i,j)*delta(i,mu)*delta(j,nu) + D2 * dAlaplacianA(i,j,mu,nu);
}

uint32_t BrusselatorBase::ip1( uint32_t i ) const {
	uint32_t r = i+1;
	if( i == nX-1 )
		r = 0;
	return r;
}

uint32_t BrusselatorBase::im1( uint32_t i ) const {
	uint32_t r = i-1;
	if( i == 0 )
		r = nX-1;
	return r;
}

uint32_t BrusselatorBase::jp1( uint32_t j ) const {
	uint32_t r = j+1;
	if( j == nY-1 )
		r = 0;
	return r;
}

uint32_t BrusselatorBase::jm1( uint32_t j ) const {
	uint32_t r = j-1;
	if( j == 0 )
		r = nY-1;
	return r;
}

double BrusselatorBase::laplacianX( const VectorXd &x, uint32_t i, uint32_t j ) const {
	return	( X(x,ip1(i),j) - 2.0*X(x,i,j) + X(x,im1(i),j) ) / ( dx*dx )
		+	( X(x,i,jp1(j)) - 2.0*X(x,i,j) + X(x,i,jm1(j)) ) / ( dy*dy );
}

double BrusselatorBase::laplacianY( const VectorXd &x, uint32_t i, uint32_t j ) const {
	return	( Y(x,ip1(i),j) - 2.0*Y(x,i,j) + Y(x,im1(i),j) ) / ( dx*dx )
		+	( Y(x,i,jp1(j)) - 2.0*Y(x,i,j) + Y(x,i,jm1(j)) ) / ( dy*dy );
}

double BrusselatorBase::dAlaplacianA( uint32_t i, uint32_t j, uint32_t mu, uint32_t nu ) const {
	return	( delta(ip1(i),mu) - 2.0*delta(i,mu) + delta(im1(i),mu) )*delta(j,nu) / ( dx*dx )
		+	( delta(jp1(j),nu) - 2.0*delta(j,nu) + delta(jm1(j),nu) )*delta(i,mu) / ( dy*dy );
}

