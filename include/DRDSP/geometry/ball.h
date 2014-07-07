#ifndef INCLUDED_GEOMETRY_BALL
#define INCLUDED_GEOMETRY_BALL
#include "../types.h"

namespace DRDSP {

	/**
	 * \brief A ball in n-dimensions.
	 */
	template<typename Scalar,int Dim>
	struct Ball {
		Matrix<Scalar,Dim,1> centre;
		Scalar radius;

		Ball() = default;

		explicit Ball( uint32_t dim ) : centre(dim), radius() {
			centre.setZero();
		}

		void Scale( Scalar factor ) {
			radius *= factor;
		}

		template<typename Derived>
		void Translate( const MatrixBase<Derived>& delta ) {
			centre += delta;
		}
		
		void SetZero() {
			centre.setZero();
			radius = Scalar();
		}

	};

	typedef Ball<double, 1> Ball1d;
	typedef Ball<double, 2> Ball2d;
	typedef Ball<double, 3> Ball3d;
	typedef Ball<double, 4> Ball4d;
	typedef Ball<double,-1> BallXd;

	template<typename Scalar>
	Scalar Volume( const Ball<Scalar,2>& ball ) {
		return M_PI * ball.radius * ball.radius;
	}

	template<typename Scalar>
	Scalar Volume( const Ball<Scalar,3>& ball ) {
		return (4.0 * M_PI / 3.0) * ball.radius * ball.radius * ball.radius;
	}

}

#endif
