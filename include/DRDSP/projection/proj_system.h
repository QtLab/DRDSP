#ifndef INCLUDED_PROJ_SYSTEM
#define INCLUDED_PROJ_SYSTEM
#include "../types.h"
#include "../data/data_set.h"
#include "proj_secant.h"
#include "../dynamics/model_orig.h"

namespace DRDSP {

	struct ProjSecantSystem {
		MatrixXd W;
		double targetMinProjectedLength,
		       minProjectedLength,
			   cost;
		uint32_t maxIterations,
		         targetDimension;

		ProjSecantSystem();
		void GetInitial();
		void Write();
		bool Read( const char* filename );
		void AnalyseSecants();
		MatrixXd Find( const SecantsSystem& secants );

		static double CostFunction( const SecantsSystem& secants, const MatrixXd &X );
		static MatrixXd CostFunctionDerivative( const SecantsSystem& secants, const MatrixXd &X );
		static double Cost( const MatrixXd &X, const void* obj );
		static MatrixXd GradCost( const MatrixXd &X, const void* obj );
	};

}

#endif
