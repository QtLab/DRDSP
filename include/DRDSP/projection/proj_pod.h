#ifndef INCLUDED_PROJ_POD
#define INCLUDED_PROJ_POD
#include <DRDSP/types.h>
#include <DRDSP/data/data_system.h>
#include <Eigen/SVD>

namespace DRDSP {

	struct ProjPOD {
		MatrixXd W;
		uint32_t targetDimension = 2;

		explicit ProjPOD( uint32_t targetDimension ) : targetDimension(targetDimension) {}
		ProjPOD& Find( const DataSet& data );
		ProjPOD& Find( const DataSystem& data );
		const ProjPOD& Write( const char* filename ) const;

	protected:
		JacobiSVD<MatrixXd>::SingularValuesType singularValues;
	};

}

#endif
