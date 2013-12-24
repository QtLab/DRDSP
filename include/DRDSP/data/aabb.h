#ifndef INCLUDED_DATA_AABB
#define INCLUDED_DATA_AABB
#include "../types.h"

namespace DRDSP {

	/*!
	 * \brief An axis-aligned bounding box in n-dimensions.
	 */
	struct AABB {
		VectorXd bMin, //!< minimum bounds
			     bMax; //!< maximum bounds

		AABB( uint32_t dim );
		double Volume() const;
		void Scale( double factor );
		void Translate( const VectorXd& delta );
		void SetZero();
	};

}

#endif
