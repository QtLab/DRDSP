#ifndef INCLUDED_PROJ_SECANT
#define INCLUDED_PROJ_SECANT
#include "../types.h"
#include "../data/data_set.h"
#include "../data/secants.h"

namespace DRDSP {

	struct SecantsSystem {
		const Secants* secants;
		uint32_t N;
	};

	struct ProjSecant {
		MatrixXd W;
		double targetMinProjectedLength;
		uint32_t maxIterations,
		         targetDimension;

		ProjSecant();
		void GetInitial( const DataSet& data );
		void GetInitial( const DataSystem& data );
		MatrixXd Find( const Secants& secants );
		MatrixXd Find( const Secants* secants, uint32_t N );
		void AnalyseSecants( const Secants* secants, uint32_t N ) const;
		void AnalyseSecants( const Secants& secants ) const;
		void Write();
		bool Read( const char* filename );

		static double CostFunction( const Secants& secants, const MatrixXd &X );
		static double CostFunction( const Secants* secants, uint32_t N, const MatrixXd &X );
		static MatrixXd CostFunctionDerivative( const Secants& secants, const MatrixXd &X );
		static MatrixXd CostFunctionDerivative( const Secants* secants, uint32_t N, const MatrixXd &X );
		static double Cost( const MatrixXd &X, const void* obj );
		static double CostN( const MatrixXd &X, const void* obj );
		static MatrixXd GradCost( const MatrixXd &X, const void* obj );
		static MatrixXd GradCostN( const MatrixXd &X, const void* obj );
	};

}

#endif
