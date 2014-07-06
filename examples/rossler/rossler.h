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

	Vector3d ComputeLinear( const VectorXd& x ) const {
		Vector3d L;
		L[0] = 0.0;
		L[1] = 0.0;
		L[2] = -x[2];
		return L;
	}

	Vector3d ComputeTranslation( const VectorXd& x ) const {
		RosslerModel rossler;
		Vector3d T;
		T[0] = -x[1] - x[2];
		T[1] = x[0] + rossler.a * x[1];
		T[2] = rossler.b + x[0] * x[2];
		return T;
	}

	Matrix<double,9,1> ComputeLinearDerivative( const VectorXd& ) const {
		Matrix<double,9,1> L;
		L.setZero();
		L[8] = -1.0;
		return L;
	}

	Matrix<double,9,1> ComputeTranslationDerivative( const VectorXd& x ) const {
		RosslerModel rossler;
		Matrix<double,9,1> T;
		T[0] = 0.0;
		T[1] = 1.0;
		T[2] = x[2];
		T[3] = -1.0;
		T[4] = rossler.a;
		T[5] = 0.0;
		T[6] = -1.0;
		T[7] = 0.0;
		T[8] = x[0];
		return T;
	}
};

struct RosslerHighModel : Model<> {
	
	RosslerHighModel() : RosslerHighModel(100) {}
	
	explicit RosslerHighModel( uint32_t n ) : Model<>(n) {}

	template<typename Derived>
	Matrix<typename Derived::Scalar,-1,1> operator()( const MatrixBase<Derived>& x ) const {
		Matrix<Derived::Scalar,-1,1> xDot;
		xDot.setZero(stateDim);
		xDot.head<3>() = rossler( (Ainv * x).head<3>() );
		return A * xDot;
	}

	MatrixXd Partials( const VectorXd& x ) const {
		return AutoDerivative( *this, x );
	}

	friend struct RosslerHighFamily;
	
	static MatrixXd A, Ainv;
	
	RosslerModel rossler;

	static void GenerateA( uint32_t dimension ) {
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
		A = svd.matrixU() * svalues.asDiagonal() * svd.matrixV().adjoint();
		Ainv = svd.matrixV() * svalues.asDiagonal().inverse() * svd.matrixU().adjoint();
	}

};

struct RosslerHighFamily : Family<RosslerHighModel> {
	
	RosslerHighFamily() : RosslerHighFamily(100) {};

	explicit RosslerHighFamily( uint32_t n ) : Family<RosslerHighModel>(n,1) {}
	
	RosslerHighModel operator()( const VectorXd& parameter ) const {
		RosslerHighModel model( stateDim );
		model.rossler.c = parameter[0];
		return model;
	}
};

#endif
