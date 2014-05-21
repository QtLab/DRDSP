#ifndef INCLUDED_PROJ_SECANT
#define INCLUDED_PROJ_SECANT
#include "../types.h"
#include "../data/data_system.h"
#include "../data/secants.h"

namespace DRDSP {

	struct SecantCostFunction {
		const SecantsPreComputed& secants;

		explicit SecantCostFunction( const SecantsPreComputed& secants ) : secants(secants) {}
		double operator()( const MatrixXd& X ) const;
	};

	struct SecantCostGradient {
		const SecantsPreComputed& secants;

		explicit SecantCostGradient( const SecantsPreComputed& secants ) : secants(secants) {}
		MatrixXd operator()( const MatrixXd& X ) const;
	};

	struct SecantCostFunctionMulti {
		const SecantsPreComputed* secants;
		size_t N;

		SecantCostFunctionMulti( const SecantsPreComputed* secants, size_t N ) : secants(secants), N(N) {}
		double operator()( const MatrixXd& X ) const;
	};

	struct SecantCostGradientMulti {
		const SecantsPreComputed* secants;
		size_t N;

		SecantCostGradientMulti( const SecantsPreComputed* secants, size_t N ) : secants(secants), N(N) {}
		MatrixXd operator()( const MatrixXd& X ) const;
	};

	struct ProjSecant {
		MatrixXd W;
		double targetMinProjectedLength;
		uint32_t targetDimension, maxIterations;

		ProjSecant();
		explicit ProjSecant( uint32_t targetDimension );
		ProjSecant& ComputeInitial( const DataSet& data );
		ProjSecant& ComputeInitial( const DataSystem& data );
		ProjSecant& Find( const SecantsPreComputed& secants );
		ProjSecant& Find( const SecantsPreComputed* secants, size_t N );
		ProjSecant& Find( const vector<SecantsPreComputed>& secants );
		const ProjSecant& AnalyseSecants( const SecantsPreComputed& secants ) const;
		const ProjSecant& AnalyseSecants( const SecantsPreComputed* secants, size_t N ) const;
		const ProjSecant& AnalyseSecants( const vector<SecantsPreComputed>& secants ) const;
		const ProjSecant& WriteCSV( const char* filename ) const;
		const ProjSecant& WriteBinary( const char* filename ) const;
		bool ReadBinary( const char* filename );
	};

}

#endif
