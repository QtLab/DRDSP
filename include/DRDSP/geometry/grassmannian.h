#ifndef INCLUDED_GEOMETRY_GRASSMANNIAN
#define INCLUDED_GEOMETRY_GRASSMANNIAN
#pragma warning( disable : 4714 ) // function '...' marked as __forceinline not inlined
#include <Eigen/SVD>
#include "metric.h"

using namespace Eigen;

namespace DRDSP {

	namespace Grassmannian {

		MatrixXd HorizontalComponent( const MatrixXd& W, const MatrixXd& V );

		MatrixXd VerticalComponent( const MatrixXd& W, const MatrixXd& V );

		using Metric = UniformMetric<MatrixXd,FrobeniusInnerProduct>;

		struct Geodesic {
			typedef Metric Metric;
			typedef double T;

			MatrixXd position, velocity;

			Geodesic() = default;

			Geodesic( const MatrixXd& x, const MatrixXd& v ) :
				position(x),
				velocity(v),
				svd( v, ComputeThinU | ComputeThinV )
			{}

			void Set( const MatrixXd& x, const MatrixXd& v );
		
			MatrixXd operator()( double t ) const;
		
			MatrixXd ParallelTranslate( const MatrixXd& V, double t ) const;

		protected:
			JacobiSVD<MatrixXd> svd;
			typedef JacobiSVD<MatrixXd>::SingularValuesType svType;
		};

	}

}

#endif
