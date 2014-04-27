#ifndef INCLUDED_GEOMETRY_METRIC
#define INCLUDED_GEOMETRY_METRIC
#include <Eigen/Core>

using namespace Eigen;

namespace DRDSP {

	template<typename V>
	struct Traits;

	template<typename V>
	struct DotProduct {
		typedef typename V Vector;
		typedef typename Traits<V>::Scalar Scalar;

		Scalar operator()( const V& x, const V& y ) const {
			Scalar r(0);
			for(uint32_t i=0;i<Traits<V>::dimension;++i) {
				r += x[i] * y[i];
			}
			return r;
		}

	};

	template<typename Point,typename InnerProduct>
	struct UniformMetric {
		typedef Point Point;
		typedef InnerProduct InnerProduct;
		InnerProduct operator()( const Point& x ) const {
			return InnerProduct();
		}
	};

	struct FrobeniusInnerProduct {
		typedef MatrixXd Vector;
		typedef MatrixXd::Scalar Scalar;

		Scalar operator()( const MatrixXd& x, const MatrixXd& y ) const {
			return ( y.adjoint() * x ).trace();
		}

	};

}

#endif
