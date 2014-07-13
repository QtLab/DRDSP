#ifndef INCLUDED_GOODFELLOW
#define INCLUDED_GOODFELLOW
#include <DRDSP/dynamics/model.h>
#include <DRDSP/auto_diff.h>

using namespace DRDSP;

struct Goodfellow : Model<> {
	MatrixXd adjacency;
	VectorXd mu;
	double omega, a, b, c, d, p;
	uint32_t N;

	Goodfellow();
	explicit Goodfellow( uint32_t numCompartments );
	void Create( uint32_t numCompartments );

	MatrixXd Partials( const VectorXd& x ) const {
		return AutoDerivative( *this, x );
	}

	template<typename Derived>
	Matrix<typename Derived::Scalar,-1,1> operator()( const MatrixBase<Derived>& state ) const {
		typedef typename Derived::Scalar Scalar;
		Matrix<Scalar,-1,1> res(2*N);
		for(uint32_t i=0;i<N;++i) {
			Scalar temp1 = omega - d * r2(state,i);
			Scalar temp2 = poly(state,i);
			res[i] = y(state,i) * temp1 + x(state,i) * temp2;
			res[N+i] = -x(state,i) * temp1 + y(state,i) * temp2;
		}
		res.head(N) += (p / N) * (adjacency * state.head(N));
		return res;
	}

	template<typename Derived>
	typename Derived::Scalar x( const MatrixBase<Derived>& state, uint32_t i ) const {
		return state[i];
	}

	template<typename Derived>
	typename Derived::Scalar y( const MatrixBase<Derived>& state, uint32_t i ) const {
		return state[N+i];
	}

	template<typename Derived>
	typename Derived::Scalar r2( const MatrixBase<Derived>& state, uint32_t i ) const {
		return x(state,i)*x(state,i) + y(state,i)*y(state,i);
	}

	template<typename Derived>
	typename Derived::Scalar r4( const MatrixBase<Derived>& state, uint32_t i ) const {
		return r2(state,i)*r2(state,i);
	}

	template<typename Derived>
	typename Derived::Scalar r6( const MatrixBase<Derived>& state, uint32_t i ) const {
		return r4(state,i)*r2(state,i);
	}

	template<typename Derived>
	typename Derived::Scalar poly( const MatrixBase<Derived>& state, uint32_t i ) const {
		return mu[i] - a * r2(state,i) + b * r4(state,i) - c * r6(state,i);
	}

};

struct GoodfellowFamily : Family<Goodfellow> {
	uint32_t N;

	explicit GoodfellowFamily( uint32_t N ) : Family<Goodfellow>(2*N,1), N(N) {}

	Goodfellow operator()( const VectorXd& parameter ) const {
		Goodfellow model( N );
		model.p = parameter[0];
		return model;
	}
};

#endif
