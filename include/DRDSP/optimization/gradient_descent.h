#ifndef INCLUDED_OPTIMIZATION_GRADIENT_DESCENT
#define INCLUDED_OPTIMIZATION_GRADIENT_DESCENT
#include "linesearch.h"

namespace DRDSP {
	
	template<typename Geodesic,typename S,typename DS>
	struct GradientDescent {
		typedef typename Geodesic::Metric Metric;
		typedef typename Metric::Point Point;

		const S& cost;
		const DS& gradient;
		Metric metric;
		LineSearch lineSearch;
		uint32_t maxSteps;
		
		GradientDescent( const S& cost, const DS& gradient ) :
			cost(cost),
			gradient(gradient),
			maxSteps(1000)
		{}
		
		bool Step( Point& x ) {

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

		bool Optimize( Point& x ) {
			uint32_t n = 0;
			if( lineSearch.alpha <= 0.0 )
				lineSearch.alpha = 1.0e-2;
			for(uint32_t i=0;i<maxSteps;++i) {
				if( !Step(x) )
					return true;
				++n;
				cout << n << "\t" << lineSearch.S0 << "\t" << lineSearch.alpha << endl;
			}
			return false;
		}

	};

}

#endif

