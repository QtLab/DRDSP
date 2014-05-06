#ifndef INCLUDED_PROJ_SECANT
#define INCLUDED_PROJ_SECANT
#include "../types.h"
#include "../data/data_system.h"
#include "../data/secants.h"

namespace DRDSP {

	struct SecantCostFunction {
		const SecantsPreComputed& secants;

		SecantCostFunction( const SecantsPreComputed& secants ) : secants(secants) {}
		double operator()( const MatrixXd& X ) const;
	};

	struct SecantCostGradient {
		const SecantsPreComputed& secants;

		SecantCostGradient( const SecantsPreComputed& secants ) : secants(secants) {}
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
		uint32_t maxIterations, targetDimension;

		ProjSecant();
		void GetInitial( const DataSet& data );
		void GetInitial( const DataSystem& data );
		void Find( const SecantsPreComputed& secants );
		void Find( const SecantsPreComputed* secants, size_t N );
		void Find( const vector<SecantsPreComputed>& secants );
		void AnalyseSecants( const SecantsPreComputed& secants ) const;
		void AnalyseSecants( const SecantsPreComputed* secants, size_t N ) const;
		void AnalyseSecants( const vector<SecantsPreComputed>& secants ) const;
		void WriteCSV( const char* filename ) const;
		void WriteBinary( const char* filename ) const;
		bool ReadBinary( const char* filename );
	};

}

#endif
