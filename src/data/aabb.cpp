#include <DRDSP/data/aabb.h>

using namespace DRDSP;

AABB::AABB( uint32_t dim ) : bMin(dim), bMax(dim) {}
		
double AABB::Volume() const {
	return (bMax - bMin).cwiseAbs().prod();
}

void AABB::Scale( double factor ) {
	VectorXd middle = (bMax + bMin) * 0.5;
	VectorXd halfDiff = (bMax - bMin) * 0.5;
	bMin = middle - halfDiff * factor;
	bMax = middle + halfDiff * factor;
}

void AABB::Translate( const VectorXd& delta ) {
	bMin += delta;
	bMax += delta;
}

void AABB::SetZero() {
	bMin.setZero();
	bMax.setZero();
}
