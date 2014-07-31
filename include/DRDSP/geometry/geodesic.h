#ifndef INCLUDED_GEOMETRY_GEODESIC
#define INCLUDED_GEOMETRY_GEODESIC

namespace DRDSP {

	template<typename V>
	struct Traits;

	template<typename V>
	struct EuclideanGeodesic {
		typedef V Point;
		typedef typename Traits<V>::Scalar T;
		
		V position, velocity;

		V operator()( T t ) const {
			return position + velocity * t;
		}

		V ParallelTranslate( const V& v, T ) const {
			return v;
		}

	};

}

#endif
