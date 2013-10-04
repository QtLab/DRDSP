#ifndef INCLUDED_OPTIMIZATION_CONJUGATEGRADIENT
#define INCLUDED_OPTIMIZATION_CONJUGATEGRADIENT
#include "optimization.h"
#include "../geometry/metric.h"

namespace DRDSP {
	template<typename TGeo,typename TMet>
	struct ConjugateGradient : Optimization<TGeo,TMet> {
		enum modType { FR, PR, HS } modifier;

		ConjugateGradient() : Optimization<TGeo,TMet>(), modifier(HS) {}

		bool Step( TVec &X ) {
			lastGrad = grad;
			lineSearch.gradSx = lineSearch.gradS(X,lineSearch.obj);
			grad = lineSearch.gradSx;

			TVec Lambda = -grad;
			if( n > 0 )
				Lambda += geodesic.ParallelTranslate(geodesic.GetTangent(),lineSearch.alpha) * Modifier(X);

			geodesic.Set(X,Lambda);

			double alpha = lineSearch.Search();

			if( alpha <= 0.0 )
				return false;

			lastX = X;
			X = geodesic.Evaluate(alpha);
			n++;

			return true;
		}
	protected:
		TVec lastX, grad, lastGrad;

		double Modifier( const TVec &X ) {
			if( n==0 ) return 0.0;
			const TMet &M = *(lineSearch.metric);
			TVec ptlastGrad, GmptlastG;
			TMet::TIP IPX = M(X);
			switch( modifier ) {
				case FR:
					return IPX.Norm2(grad) / M(lastX).Norm2(lastGrad);
				break;
				case PR:
					ptlastGrad = geodesic.ParallelTranslate(lastGrad,lineSearch.alpha);
					return  IPX(grad,grad-ptlastGrad) / IPX.Norm2(ptlastGrad);
				break;
				case HS:
					ptlastGrad = geodesic.ParallelTranslate(lastGrad,lineSearch.alpha);
					GmptlastG = grad-ptlastGrad;
					return IPX(grad,GmptlastG) / IPX(ptlastGrad,GmptlastG);
				break;
			}
			return 0.0;
		}
	};
}

#endif

