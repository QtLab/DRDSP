#ifndef INCLUDED_PROJ_POD
#define INCLUDED_PROJ_POD
#include <DRDSP/types.h>
#include <DRDSP/data/data_system.h>
#include <Eigen/SVD>

namespace DRDSP {

	struct ProjPOD {
		MatrixXd W;
		uint32_t targetDimension;

		ProjPOD();
		explicit ProjPOD( uint32_t targetDimension );
		void Find( const DataSet& data );
		void Find( const DataSystem& data );
		void Write( const char* filename ) const;

	protected:
		MatrixXd dataMatrix;
		JacobiSVD<MatrixXd> svd;
	};

}

#endif
