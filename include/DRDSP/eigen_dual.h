#ifndef INCLUDED_EIGEN_DUAL
#define INCLUDED_EIGEN_DUAL
#include "dual.h"
#include <Eigen/Core>

using namespace Eigen;

namespace DRDSP {

	template<typename T, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
	using DualMatrix = Matrix<dual<T>,Rows,Cols,Options,MaxRows,MaxCols>;
	
	typedef DualMatrix<float,-1,-1,0,-1,-1>  DualMatrixXf;
	typedef DualMatrix<double,-1,-1,0,-1,-1> DualMatrixXd;

}

#endif
