#ifndef INCLUDED_DATA_AABB
#define INCLUDED_DATA_AABB
#include "../types.h"

namespace DRDSP {

	struct AABB {
		VectorXd bMin, bMax;

		AABB( uint32_t dim );
		double Volume() const;
		void Scale( double factor );
		void Translate( const VectorXd& delta );
		void SetZero();
	};

}

#endif
