#ifndef INCLUDED_OPTIMIZATION_OPTIMIZATION
#define INCLUDED_OPTIMIZATION_OPTIMIZATION
#include <iostream>
#include "../types.h"
#include "linesearch.h"
#include "../Geometry/geodesic.h"

using namespace std;

namespace DRDSP {
	template<typename TGeo,typename TMet>
	struct Optimization {
		typedef typename TMet::TVec TVec;
		LineSearch<TMet> lineSearch;
		uint32_t n, maxSteps;

		Optimization() : lineSearch(&geodesic), n(0), maxSteps(0) {}

		virtual bool Step( TVec &X ) = 0;

		virtual bool Optimize( TVec &X ) {
			n = 0;
			if( lineSearch.alpha <= 0.0 ) lineSearch.alpha = 1.0e-2;
			for(uint32_t i=0;i<maxSteps;i++) {
				if( !Step(X) ) return true;
				cout << n << "\t" << lineSearch.Sx << "\t" << lineSearch.alpha << endl;
			}
			return false;
		}

	protected:
		TGeo geodesic;

	};
}

#endif

