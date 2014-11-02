#ifndef INCLUDED_GEOMETRY_GRASSMANNIAN
#define INCLUDED_GEOMETRY_GRASSMANNIAN
#include "metric.h"

#pragma warning( disable : 4510 ) // default constructor could not be generated
#pragma warning( disable : 4610 ) // can never be instantiated - user defined constructor required

#include <Eigen/SVD>

#pragma warning( default : 4610 )
#pragma warning( default : 4510 )

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
