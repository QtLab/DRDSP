#include <DRDSP/dynamics/embedding.h>

using namespace DRDSP;

Embedding::Embedding( uint32_t origDim, uint32_t embedDim ) : oDim(origDim), eDim(embedDim) {}

EmbeddingCW::EmbeddingCW( uint32_t origDim, uint32_t embedDim ) : oDim(origDim), eDim(embedDim) {}

VectorXd Embedding::Evaluate( const VectorXd &x ) const {
	return x;
}

VectorXd EmbeddingCW::Evaluate( const VectorXd &x ) const {
	return x;
}

MatrixXd Embedding::Derivative( const VectorXd &x ) const {
	return MatrixXd::Identity(eDim,oDim);
}

double EmbeddingCW::Derivative( const VectorXd &x, uint32_t i, uint32_t j ) const {
	return (i==j)?1.0:0.0;
}

MatrixXd Embedding::DerivativeAdjoint( const VectorXd &x ) const {
	return MatrixXd::Identity(oDim,eDim);
}

double EmbeddingCW::DerivativeAdjoint( const VectorXd &x, uint32_t i, uint32_t j ) const {
	return (i==j)?1.0:0.0;
}

MatrixXd Embedding::Derivative2( const VectorXd &x, uint32_t i ) const {
	return MatrixXd::Identity(oDim,oDim);
}

double EmbeddingCW::Derivative2( const VectorXd &x, uint32_t i, uint32_t j, uint32_t k ) const {
	return (j==k)?1.0:0.0;
}

MatrixXd Embedding::ComputeInducedMetric( const VectorXd &x ) const {
	MatrixXd deriv = Derivative(x);
	return deriv.transpose() * deriv;
}

double EmbeddingCW::ComputeInducedMetric( const VectorXd &x, uint32_t i, uint32_t j ) const {
	double r = 0.0;
	if( i==j ) {
		double temp;
		for(uint32_t a=0;a<eDim;a++) {
			temp = Derivative(x,a,i);
			r += temp * temp;
		}
	} else {
		for(uint32_t a=0;a<eDim;a++)
			r += Derivative(x,a,i) * Derivative(x,a,j);
	}
	return r;
}

DataSet Embedding::EmbedData( const DataSet& data ) const {
	DataSet r( data.points.size(), data.dimension );
	for(uint32_t i=0;i<r.points.size();i++) {
		r[i] = Evaluate(data[i]);
	}
	return r;
}

DataSystem Embedding::EmbedData( const DataSystem& data ) const {
	DataSystem r( eDim, data.numParameters, data.parameterDimension );
	for(uint32_t i=0;i<r.numParameters;i++) {
		r.parameters[i] = data.parameters[i];
		r.dataSets[i].points.resize( data.dataSets[i].points.size() );

		for(uint32_t j=0;j<r.dataSets[i].points.size();j++)
			r.dataSets[i][j] = Evaluate(data.dataSets[i][j]);
	}

	return r;
}

DataSet EmbeddingCW::EmbedData( const DataSet& data ) const {
	DataSet r( data.points.size(), data.dimension );
	for(uint32_t i=0;i<r.points.size();i++) {
		r[i] = Evaluate(data[i]);
	}
	return r;
}

DataSystem EmbeddingCW::EmbedData( const DataSystem& data ) const {
	DataSystem r( data.dimension, data.numParameters, data.parameterDimension );
	for(uint32_t i=0;i<r.numParameters;i++)
		for(uint32_t j=0;j<r.dataSets[i].points.size();j++) {
			r.dataSets[i][j] = Evaluate(data.dataSets[i][j]);
		}
	return r;
}