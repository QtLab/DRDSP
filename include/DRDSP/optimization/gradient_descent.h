#ifndef INCLUDED_OPTIMIZATION_GRADIENTDESCENT
#define INCLUDED_OPTIMIZATION_GRADIENTDESCENT
#include "optimization.h"

namespace DRDSP {
	template<typename TGeo,typename TMet>
	struct GradientDescent : Optimization<TGeo,TMet> {
		bool Step( TVec &X ) {
			lineSearch.gradSx = lineSearch.gradS(X,lineSearch.obj);
			geodesic.Set(X,-lineSearch.gradSx);

			double alpha = lineSearch.Search();

			if( alpha <= 0.0 )
				return false;

			X = geodesic.Evaluate(alpha);
			n++;
			return true;
		}
	};
}

#endif

