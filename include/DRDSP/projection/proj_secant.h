#ifndef INCLUDED_PROJ_SECANT
#define INCLUDED_PROJ_SECANT
#include "../types.h"
#include "../data/data_system.h"
#include "../data/secants.h"

namespace DRDSP {

	struct SecantCostFunction {
		const Secants& secants;

		explicit SecantCostFunction( const Secants& secants ) : secants(secants) {}
		double operator()( const MatrixXd& X ) const;
	};

	struct SecantCostGradient {
		const Secants& secants;

		explicit SecantCostGradient( const Secants& secants ) : secants(secants) {}
		MatrixXd operator()( const MatrixXd& X ) const;
	};

	struct SecantCostFunctionMulti {
		const Secants* secants;
		size_t N;

		SecantCostFunctionMulti( const Secants* secants, size_t N ) : secants(secants), N(N) {}
		double operator()( const MatrixXd& X ) const;
	};

	struct SecantCostGradientMulti {
		const Secants* secants;
		size_t N;

		SecantCostGradientMulti( const Secants* secants, size_t N ) : secants(secants), N(N) {}
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
		ProjSecant& Find( const Secants& secants );
		ProjSecant& Find( const Secants* secants, size_t N );
		ProjSecant& Find( const vector<Secants>& secants );
		const ProjSecant& AnalyseSecants( const Secants& secants ) const;
		const ProjSecant& AnalyseSecants( const Secants* secants, size_t N ) const;
		const ProjSecant& AnalyseSecants( const vector<Secants>& secants ) const;
		const ProjSecant& WriteCSV( const char* filename ) const;
		const ProjSecant& WriteBinary( const char* filename ) const;
		bool ReadBinary( const char* filename );
	};

}

#endif
