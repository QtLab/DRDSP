#ifndef INCLUDED_GEOMETRY_INNER_PRODUCT
#define INCLUDED_GEOMETRY_INNER_PRODUCT
#include "../types.h"
#include "vector.h"

namespace DRDSP {	
	template<typename Derived>
	struct InnerProduct {
		typedef typename Traits<Derived>::TVec TVec;
		typedef typename Traits<Derived>::TCV TCV;
		typedef typename Traits<Derived>::TScalar TScalar;
		virtual TScalar operator()( const TVec& V, const TVec& W ) const = 0;
		virtual TScalar Norm2( const TVec& V ) const = 0;
		virtual TCV Flat( const TVec& V ) const = 0;
		virtual TVec Sharp( const TCV& V ) const = 0;
	};

	template<>
	struct Traits<struct IPFrobenius> {
		typedef CoMatrixXd TCV;
		typedef Traits<TCV>::TVec TVec;
		typedef Traits<TVec>::TScalar TScalar;
	};

	struct IPFrobenius : InnerProduct<IPFrobenius> {
		TScalar operator()( const TVec& V, const TVec& W ) const {
			return ( W.adjoint() * V ).trace();
		}
		TScalar Norm2( const TVec& V ) const {
			return operator()(V,V);
		}
		TCV Flat( const TVec& V ) const {
			return V.adjoint();
		}
		TVec Sharp( const TCV& V ) const {
			return V.adjoint();
		}
	};

	template<>
	struct Traits<struct IPVector> {
		typedef CoVectorXd TCV;
		typedef Traits<TCV>::TVec TVec;
		typedef Traits<TVec>::TScalar TScalar;
	};

	template<>
	struct Traits<struct IPDot> {
		typedef CoVectorXd TCV;
		typedef Traits<TCV>::TVec TVec;
		typedef Traits<TVec>::TScalar TScalar;
	};

	struct IPDot : InnerProduct<IPDot> {
		TScalar operator()( const TVec& V, const TVec& W ) const {
			return V.dot( W );
		}
		TScalar Norm2( const TVec& V ) const {
			return operator()(V,V);
		}
		TCV Flat( const TVec& V ) const {
			return V;
		}
		TVec Sharp( const TCV& V ) const {
			return V;
		}
	};
}

#endif

