#ifndef INCLUDED_ROSSLER
#define INCLUDED_ROSSLER
#include <DRDSP/dynamics/model.h>
#include <DRDSP/auto_diff.h>
#include <Eigen/SVD>

using namespace DRDSP;

struct RosslerModel : Model<Vector3d> {
	double a,b,c;

	RosslerModel() : RosslerModel(4) {}
	
	explicit RosslerModel( double c ) : Model<Vector3d>(3), a(0.1), b(0.1), c(c) {}

	template<typename Derived>
	Matrix<typename Derived::Scalar,3,1> operator()( const MatrixBase<Derived>& x ) const {
		Matrix<typename Derived::Scalar,3,1> v;
		v[0] = -x[1] - x[2];
		v[1] = x[0] + a * x[1];
		v[2] = b + x[2] * ( x[0] - c );
		return v;
	}

	Matrix3d Partials( const Vector3d& x ) const {
		return AutoDerivative( *this, x );
	}
};

struct RosslerFamily : Family<RosslerModel> {
	
	RosslerFamily() : Family<RosslerModel>(3,1) {}
	
	RosslerModel operator()( const VectorXd& parameter ) const {
		return RosslerModel( parameter[0] );
	}

};

struct RosslerHighModel : Model<> {
	
	RosslerHighModel() : RosslerHighModel(100) {}
	
	explicit RosslerHighModel( uint32_t n ) : Model<>(n) {
		GenerateA();
	}

	template<typename Derived>
	Matrix<typename Derived::Scalar,-1,1> operator()( const MatrixBase<Derived>& x ) const {
		Matrix<Derived::Scalar,-1,1> xDot;
		xDot.setZero(dimension);
		xDot.head<3>() = model( (Ainv * x).head<3>() );
		return A * xDot;
	}

	Matrix3d Partials( const Vector3d& x ) const {
		return AutoDerivative( *this, x );
	}

protected:
	RosslerModel model;
	MatrixXd A, Ainv;

	void GenerateA() {
		MatrixXd R(dimension,dimension);
		R.setRandom();
		R *= 2.0;
		R.array() -= 1.0;
		JacobiSVD<MatrixXd,NoQRPreconditioner> svd(R,ComputeFullU|ComputeFullV);
		VectorXd svalues(dimension);
		double smin = 0.5, smax = 1.5;
		for(uint32_t i=0;i<dimension;++i) {
			svalues[i] = smax + (( smin - smax )/(dimension-1)) * i;
		}
		A = svd.matrixU() * svalues.asDiagonal() * svd.matrixV().transpose();
		Ainv = svd.matrixV() * svalues.asDiagonal().inverse() * svd.matrixU().transpose();
	}

};

#endif
