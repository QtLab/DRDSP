#ifndef INCLUDED_OPTIMIZATION_GRADIENT_DESCENT
#define INCLUDED_OPTIMIZATION_GRADIENT_DESCENT
#include "linesearch.h"

namespace DRDSP {
	
	template<typename Geodesic>
	struct GradientDescent {
		typedef typename Geodesic::Metric Metric;
		typedef typename Metric::Point Point;

		Metric metric;
		LineSearch lineSearch;
		uint32_t maxSteps = 1000;
		
		template<typename S,typename DS>
		bool Step( Point& x, const S& cost, const DS& gradient ) {

			Geodesic geodesic( x, -gradient(x) );

			double alpha = lineSearch.Search(
				[&]( double t ){ return cost( geodesic(t) ); },
				[&]( double t ){
					Point xt = geodesic(t);
					return metric(xt)(
						gradient(xt),
						geodesic.ParallelTranslate( geodesic.velocity, t )
					);
				}
			);

			if( alpha <= 0.0 ) return false;

			x = geodesic(alpha);

			return true;
		}

		template<typename S,typename DS>
		bool Optimize( Point& x, const S& cost, const DS& gradient ) {
			uint32_t n = 0;
			if( lineSearch.alpha <= 0.0 )
				lineSearch.alpha = 1.0e-2;
			for(uint32_t i=0;i<maxSteps;++i) {
				if( !Step( x, cost, gradient ) )
					return true;
				++n;
				cout << n << "\t" << lineSearch.S0 << "\t" << lineSearch.alpha << endl;
			}
			return false;
		}

	};
}

#endif

