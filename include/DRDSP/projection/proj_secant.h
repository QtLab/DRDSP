#ifndef INCLUDED_PROJ_SECANT
#define INCLUDED_PROJ_SECANT
#include "../types.h"
#include "../data/data_set.h"
#include "../data/secants.h"

namespace DRDSP {

	struct SecantsSystem {
		const Secants* secants;
		uint16_t N;
	};

	struct ProjSecant {
		MatrixXd W;
		double targetMinProjectedLength;
		uint16_t maxIterations, targetDimension;

		ProjSecant();
		void GetInitial( const DataSet& data );
		void GetInitial( const DataSystem& data );
		void Find( const Secants& secants );
		void Find( const Secants* secants, uint16_t N );
		void AnalyseSecants( const Secants* secants, uint16_t N ) const;
		void AnalyseSecants( const Secants& secants ) const;
		void WriteText( const char* filename ) const;
		void WriteBinary( const char* filename ) const;
		bool ReadBinary( const char* filename );

		static double CostFunction( const Secants& secants, const MatrixXd &X );
		static double CostFunction( const Secants* secants, uint16_t N, const MatrixXd &X );
		static MatrixXd CostFunctionDerivative( const Secants& secants, const MatrixXd &X );
		static MatrixXd CostFunctionDerivative( const Secants* secants, uint16_t N, const MatrixXd &X );
		static double Cost( const MatrixXd &X, const void* obj );
		static double CostN( const MatrixXd &X, const void* obj );
		static MatrixXd GradCost( const MatrixXd &X, const void* obj );
		static MatrixXd GradCostN( const MatrixXd &X, const void* obj );
	};

}

#endif
