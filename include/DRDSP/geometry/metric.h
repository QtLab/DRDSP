#ifndef INCLUDED_GEOMETRY_METRIC
#define INCLUDED_GEOMETRY_METRIC
#include "../types.h"
#include "inner_product.h"

namespace DRDSP {
	template<typename Derived>
	struct Metric {
		typedef typename Traits<Derived>::TIP TIP;
		typedef typename Traits<TIP>::TVec TVec;
		virtual TIP operator()( const TVec& X ) const  = 0;
	};

	template<typename _TIP> struct MetricUniform;

	template<typename _TIP>
	struct Traits<MetricUniform<_TIP>> {
		typedef _TIP TIP;
		typedef typename Traits<TIP>::TCV TCV;
		typedef typename Traits<TIP>::TVec TVec;
		typedef typename Traits<TIP>::TScalar TScalar;
	};

	template<typename _TIP>
	struct MetricUniform : Metric<MetricUniform<_TIP>> {
		TIP IP;
		TIP operator()( const TVec& X ) const  {
			return IP;
		}
		const TIP& operator()() const {
			return IP;
		}
	};

	typedef MetricUniform<IPFrobenius> MetricFrobenius;
	typedef MetricUniform<IPDot> MetricDot;
}

#endif

