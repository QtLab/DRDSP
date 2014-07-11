#include <DRDSP/projection/inverse.h>
#include <DRDSP/misc.h>
#include <Eigen/LU>
#include <iostream>

using namespace std;
using namespace DRDSP;

static VectorXd MeanVertical( const vector<VectorXd>& points, const MatrixXd& W ) {
	VectorXd mean = accumulate( points, VectorXd(VectorXd::Zero(W.rows())) ) / (double)points.size();
	return mean - W * ( W.adjoint() * mean );
}

static vector<VectorXd> MeanVerticals( const DataSystem& data, const MatrixXd& W ) {
	vector<VectorXd> verticals( data.numParameters );
	for(uint32_t i=0;i<data.numParameters;++i) {
		verticals[i] = MeanVertical( data.dataSets[i].points, W );
	}
	return verticals;
}

static MatrixXd ComputeLHS( const vector<VectorXd>& parameters, const VectorXd& meanParam ) {
	MatrixXd LHS;
	int64_t pdim = meanParam.size();
	LHS.setZero( parameters.size() * pdim, pdim );
	for(size_t i=0;i<parameters.size();++i) {
		auto diff = parameters[i] - meanParam;
		LHS.block(i*pdim,0,pdim,pdim) = diff * diff.transpose();
	}
	return LHS;
}

static MatrixXd ComputeRHS( const vector<VectorXd>& parameters, const vector<VectorXd>& meanVerticals, const VectorXd& meanParam ) {
	MatrixXd RHS;
	int64_t dim = meanVerticals[0].size();
	int64_t pdim = meanParam.size();
	RHS.setZero( parameters.size() * pdim, dim );
	for(size_t i=0;i<parameters.size();++i) {
		RHS.block(i*pdim,0,pdim,dim) = (parameters[i] - meanParam) * meanVerticals[i].transpose();
	}
	return RHS;
}

ProjInverse DRDSP::ComputeInverse( const MatrixXd& W, const DataSystem& data ) {
	ProjInverse psi;
	psi.Compute( W, data );
	return psi;
}

void ProjInverse::Compute( const MatrixXd& W, const DataSystem& data ) {
	vector<VectorXd> meanVerticals = MeanVerticals( data, W );
	VectorXd meanParam = accumulate( data.parameters, VectorXd(VectorXd::Zero(data.paramDim)) ) / data.numParameters;
	MatrixXd LHS = ComputeLHS( data.parameters, meanParam );
	MatrixXd RHS = ComputeRHS( data.parameters, meanVerticals, meanParam );
	this->W = W;
	vertical.sourceDim = data.paramDim;
	vertical.targetDim = data.dimension;
	vertical.linear = FullPivLU<MatrixXd>(LHS).solve(RHS).transpose();
	vertical.translation = accumulate( meanVerticals, VectorXd(VectorXd::Zero(data.dimension)) ) / data.numParameters - vertical.linear * meanParam;
}

double DRDSP::ComputeInverseCost( const ProjInverse& inverse, const DataSystem& data ) {
	double sum = 0.0;
	for(size_t i=0;i<data.numParameters;++i) {
		double innerSum = 0.0;
		for( const auto& x : data.dataSets[i].points ) {
			innerSum += ( inverse( inverse.W.adjoint() * x, data.parameters[i] ) - x ).squaredNorm();
		}
		sum += innerSum / (double)data.dataSets[i].points.size();
	}
	return sum / data.numParameters;
}

