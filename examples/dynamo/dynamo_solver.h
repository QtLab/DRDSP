#ifndef INCLUDED_DYNAMO_SOLVER
#define INCLUDED_DYNAMO_SOLVER
#include <DRDSP/types.h>
#include "dynamo.h"

namespace DRDSP {

	struct DynamoSolver : Dynamo {
		typedef double Time;
		typedef VectorXd State;
		State state;
		double dtMax;

		explicit DynamoSolver( const Dynamo& dynamo );

		void Advance( double dt );

	protected:

		MatrixXd pZero, tZero;

		MatrixXd pMinus, pPlus, pCoef1, pCoef2, pCoef3, pCoef4, pCoef5,
		         tMinus, tPlus, tCoef1, tCoef2, tCoef3, tCoef4, tCoef5,
		         tpCof1, tpCof2, tpCof3, tpCof4, tpCof5, ptCof5,
				 diffA, under, alpha;
		void Init();
		void Step( double dt );
		void npCoff( double dTime );
		void ntCoff( double dTime );
		void npStep( double dTime );
		void ntStep( double dTime );
		double F( const MatrixXd& a, uint32_t i, uint32_t j ) const;
		double Beta( const MatrixXd& a, uint32_t i, uint32_t j ) const;
		double Btheta( const MatrixXd& a, uint32_t i, uint32_t j ) const;
		double NormB2( const MatrixXd& a, const MatrixXd& b, uint32_t i, uint32_t j ) const;
		double Alpha( const MatrixXd& a, const MatrixXd& b, uint32_t i, uint32_t j ) const;
		double dsp( const MatrixXd& p, uint32_t i, uint32_t j ) const;
		double dtp( const MatrixXd& p, uint32_t i, uint32_t j ) const;
		double eta( uint32_t i ) const;

		inline double theta( uint32_t j ) const {
			return dth * j;
		}

		inline double s( uint32_t i ) const {
			return ds * i;
		}
	};

}

#endif
