#ifndef INCLUDED_GEOMETRY_GEODESIC
#define INCLUDED_GEOMETRY_GEODESIC
#include "../types.h"

namespace DRDSP {
	template<typename TVec>
	struct Geodesic {
		virtual void Set( const TVec& point, const TVec& tangent ) {
			Point = point;
			Tangent = tangent;
		}
		const TVec& GetPoint() const {
			return Point;
		}
		const TVec& GetTangent() const {
			return Tangent;
		}
		virtual TVec Evaluate( double t ) {
			if( t == 0.0 ) return Point;
			return Point + Tangent * t;
		}
		virtual TVec ParallelTranslate( const TVec &V, double t ) {
			return V;
		}

	protected:
		TVec Point;
		TVec Tangent;
	};
}

#endif

