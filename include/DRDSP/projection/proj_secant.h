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
		SecantCostFunction& operator=( const SecantCostFunction& ) = delete;
	};

	struct SecantCostGradient {
		const Secants& secants;

		explicit SecantCostGradient( const Secants& secants ) : secants(secants) {}
		MatrixXd operator()( const MatrixXd& X ) const;
		SecantCostGradient& operator=( const SecantCostGradient& ) = delete;
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
		double targetMinProjectedLength = 0.5;
		uint32_t targetDimension = 2,
		         maxIterations = 100;

		explicit ProjSecant( uint32_t targetDimension ) : targetDimension(targetDimension) {}
		ProjSecant& ComputeInitial( const DataSet& data );
		ProjSecant& ComputeInitial( const DataSystem& data );
		ProjSecant& Find( const Secants& secants );
		ProjSecant& Find( const Secants* secants, size_t N );
		ProjSecant& Find( const vector<Secants>& secants );
		const ProjSecant& AnalyseSecants( const Secants& secants ) const;
		const ProjSecant& AnalyseSecants( const vector<Secants>& secants ) const;
		const ProjSecant& WriteCSV( const char* filename ) const;
		const ProjSecant& WriteBinary( const char* filename ) const;
		bool ReadBinary( const char* filename );
	};

}

#endif
