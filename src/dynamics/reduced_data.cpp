#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <DRDSP/dynamics/reduced_data.h>

using namespace std;
using namespace Eigen;
using namespace DRDSP;

static double eps( double x ) {
	int64_t y = *((int64_t*)&x);
	y++;
	return *((double*)&y) - x;
}

ReducedData::ReducedData() : dimension(0), count(0) {}

ReducedData::ReducedData( uint16_t dim, size_t numPoints ) : dimension(0), count(0) {
	Create(dim,numPoints);
}

void ReducedData::Create( uint16_t dim, size_t numPoints ) {
	points.resize(numPoints);
	vectors.resize(numPoints);
	derivatives.resize(numPoints);
	count = numPoints;
	dimension = dim;
	for(size_t i=0;i<count;i++) {
		points[i].setZero(dimension);
		vectors[i].setZero(dimension);
		derivatives[i].setZero(dimension,dimension);
	}
}

void ReducedData::ComputeData( Model& original, const DataSet& data, const MatrixXd& W ) {
	Create( (uint16_t)W.cols(), data.points.size() );

	for(uint32_t i=0;i<count;i++) {
		//original.PrepareOptimizations( data.points[i] );
		points[i] = W.adjoint() * data.points[i];
		vectors[i] = W.adjoint() * original.VectorField( data.points[i] );
		derivatives[i] = W.adjoint() * original.Partials( data.points[i] ) * W;
	}
	scales[0] = ComputeVectorScale();
	scales[1] = ComputeDerivativeScale();
}

void ReducedData::ComputeData( ModelParameterized& original, const VectorXd& parameter, const DataSet& data, const MatrixXd& W ) {
	Create( (uint16_t)W.cols(), data.points.size() );

	for(uint32_t i=0;i<count;i++) {
		//original.PrepareOptimizations( data.points[i], parameter );
		points[i] = W.adjoint() * data.points[i];
		vectors[i] = W.adjoint() * original.VectorField( data.points[i], parameter );
		derivatives[i] = W.adjoint() * original.Partials( data.points[i], parameter ) * W;
	}
	scales[0] = ComputeVectorScale();
	scales[1] = ComputeDerivativeScale();
}

void ReducedData::ComputeData( ModelParameterizedEmbedded& original, const VectorXd& parameter, const DataSet& data, const MatrixXd& W ) {
	Create( (uint16_t)W.cols(), data.points.size() );

	static double stabilityFactor = 1.0;

	MatrixXd A;

	for(uint32_t i=0;i<count;i++) {
		//original.PrepareOptimizations( data.points[i], parameter );

		points[i] = W.adjoint() * original.embedding.Evaluate(data.points[i]);
		vectors[i] = W.adjoint() * original.VectorField( data.points[i], parameter );

		JacobiSVD<MatrixXd> svd(W.adjoint() * original.embedding.Derivative(data.points[i]),ComputeThinU);
		uint32_t rank = 0;
		double tolerance = original.model.dimension * eps(svd.singularValues()(0));
		for(int j=0;j<svd.nonzeroSingularValues();j++) {
			if( svd.singularValues()(j) > tolerance )
				rank++;
			else break;
		}
		A = svd.matrixU().leftCols(rank) * svd.matrixU().leftCols(rank).transpose();
		derivatives[i] = W.adjoint() * original.Partials( data.points[i], parameter ) * original.embedding.DerivativeAdjoint(data.points[i]) * W * A - stabilityFactor * ( MatrixXd::Identity(A.rows(),A.cols()) - A );
	}
	scales[0] = ComputeVectorScale();
	scales[1] = ComputeDerivativeScale();
}

void ReducedData::ComputeData( ModelCW& original, const DataSet& data, const MatrixXd& W ) {
	Create( (uint16_t)W.cols(), data.points.size() );

	for(uint32_t i=0;i<count;i++) {
		//original.PrepareOptimizations( data.points[i] );

		points[i] = W.adjoint() * data.points[i];

		for(uint32_t j=0;j<dimension;j++)
			for(uint32_t k=0;k<original.dimension;k++)
				vectors[i](j) = W(k,j) * original.VectorField( data.points[i], k );

		for(uint32_t c1=0;c1<original.dimension;c1++)
			for(uint32_t c2=0;c2<original.dimension;c2++) {
				double temp = original.Partials( data.points[i], c1, c2 );
				if( temp != 0.0 ) {
					for(uint16_t k=0;k<dimension;k++)
						for(uint16_t j=0;j<dimension;j++)
							derivatives[i](j,k) += W(c1,j) * temp * W(c2,k);
				}
			}
	}
	scales[0] = ComputeVectorScale();
	scales[1] = ComputeDerivativeScale();
}

void ReducedData::ComputeData( ModelParameterizedCW& original, const VectorXd& parameter, const DataSet& data, const MatrixXd& W ) {
	Create( (uint16_t)W.cols(), data.points.size() );

	for(uint32_t i=0;i<count;i++) {
		//original.PrepareOptimizations( data.points[i], parameter );

		points[i] = W.adjoint() * data.points[i];

		for(uint32_t j=0;j<dimension;j++)
			for(uint32_t k=0;k<original.dimension;k++)
				vectors[i](j) = W(k,j) * original.VectorField( data.points[i], parameter, k );

		for(uint32_t c1=0;c1<original.dimension;c1++)
			for(uint32_t c2=0;c2<original.dimension;c2++) {
				double temp = original.Partials( data.points[i], parameter, c1, c2 );
				if( temp != 0.0 ) {
					for(uint16_t k=0;k<dimension;k++)
						for(uint16_t j=0;j<dimension;j++)
							derivatives[i](j,k) += W(c1,j) * temp * W(c2,k);
				}
			}
	}
	scales[0] = ComputeVectorScale();
	scales[1] = ComputeDerivativeScale();
}

void ReducedData::ComputeData( ModelParameterizedEmbeddedCW& original, const VectorXd& parameter, const DataSet& data, const MatrixXd& W ) {
	Create( (uint16_t)W.cols(), data.points.size() );

	static double stabilityFactor = 1.0;

	MatrixXd A, dPhi, tempMatrix;
	
	for(uint32_t i=0;i<count;i++) {
		//original.PrepareOptimizations( data.points[i], parameter );
		dPhi.setZero(dimension,original.embedding.eDim);
		tempMatrix.setZero(dimension,dimension);

		points[i] = W.adjoint() * original.embedding.Evaluate(data.points[i]);
		for(uint32_t j=0;j<dimension;j++)
			for(uint32_t k=0;k<original.embedding.eDim;k++) {
				vectors[i](j) = W(k,j) * original.VectorField( data.points[i], parameter, k );
				for(uint32_t c=0;c<original.embedding.oDim;c++)
					dPhi(j,c) += W(k,j) * original.embedding.Derivative(data.points[i],k,c);
			}

		JacobiSVD<MatrixXd> svd(dPhi,ComputeThinU);
		uint32_t rank = 0;
		double tolerance = original.model.dimension * eps(svd.singularValues()(0));
		for(int j=0;j<svd.nonzeroSingularValues();j++) {
			if( svd.singularValues()(j) > tolerance )
				rank++;
			else break;
		}
		A = svd.matrixU().leftCols(rank) * svd.matrixU().leftCols(rank).transpose();
		
		for(uint32_t c1=0;c1<original.model.dimension;c1++)
			for(uint32_t c2=0;c2<original.model.dimension;c2++) {
				double temp1 = original.Partials( data.points[i], parameter, c1, c2 );
				if( temp1 == 0.0 ) continue;
				for(uint32_t c3=0;c3<original.embedding.eDim;c3++) {
					double temp2 = original.embedding.DerivativeAdjoint(data.points[i],c2,c3);
					if( temp2 == 0.0 ) continue;
					double temp = temp1 * temp2;
					for(uint16_t k=0;k<dimension;k++)
						for(uint16_t j=0;j<dimension;j++)
							tempMatrix(j,k) += W(c1,j) * temp * W(c2,k);
				}
			}
		derivatives[i] = tempMatrix * A - stabilityFactor * ( MatrixXd::Identity(A.rows(),A.cols()) - A );
	}
	scales[0] = ComputeVectorScale();
	scales[1] = ComputeDerivativeScale();
}

AABB ReducedData::ComputeBoundingBox() const {
	AABB box(dimension);
	box.SetZero();
	for(uint32_t k=0;k<count;k++) {
		for(uint16_t j=0;j<dimension;j++) {
			if( points[k](j) > box.bMax(j) )
				box.bMax(j) = points[k](j);
			if( points[k](j) < box.bMin(j) )
				box.bMin(j) = points[k](j);
		}
	}
	return std::move(box);
}

double ReducedData::ComputeVectorScale() {
	double S1 = 0.0;
	for(uint32_t i=0;i<count;i++) {
		S1 += vectors[i].squaredNorm();
	}
	return S1 / count;
}

double ReducedData::ComputeDerivativeScale() {
	double S2 = 0.0;
	for(uint32_t i=0;i<count;i++) {
		S2 += derivatives[i].squaredNorm();
	}
	return S2 / count;
}

bool ReducedData::ReadData( const char* filename ) {
	ifstream in(filename,ios::binary);
	if( !in ) {
		cout << "ReducedData::ReadData : file error" << endl;
		return false;
	}

	in.seekg(0, ios::end);
	if( (size_t)in.tellg() < sizeof(double) * (size_t)count * (size_t)dimension * (size_t)( dimension + 2 ) ) {
		cout << "ReducedData::ReadData : insufficient data in file " << filename << endl;
		return false;
	}
	in.seekg(0, ios::beg);

	for(uint32_t k=0;k<count;k++) {
		points[k].setZero(dimension);
		vectors[k].setZero(dimension);
		derivatives[k].setZero(dimension,dimension);
		in.read( (char*)&points[k](0), sizeof(double) * dimension);
		in.read( (char*)&vectors[k](0), sizeof(double) * dimension);
		in.read( (char*)&derivatives[k](0,0), sizeof(double) * dimension * dimension );
	}
	return true;
}

void ReducedData::WriteData( const char* filename ) const {
	ofstream out(filename,ios::binary);
	if( !out ) {
		cout << "ReducedData::WriteData : file error" << endl;
		return;
	}
	for(uint32_t k=0;k<count;k++) {
		out.write( (const char*)&points[k](0), sizeof(double) * dimension );
		out.write( (const char*)&vectors[k](0), sizeof(double) * dimension );
		out.write( (const char*)&derivatives[k](0,0), sizeof(double) * dimension * dimension );
	}
}

void ReducedData::WritePointsCSV( const char* filename ) const {
	ofstream out(filename);
	if( !out ) {
		cout << "ReducedData::WritePointsText : file error" << endl;
		return;
	}
	for(uint32_t i=0;i<count;i++) {
		for(uint16_t j=0;j<dimension;j++)
			out << points[i](j) << ",";
		out << endl;
	}
}

void ReducedData::WriteVectorsCSV( const char* filename ) const {
	ofstream out(filename);
	if( !out ) {
		cout << "ReducedData::WriteVectorsText : file error" << endl;
		return;
	}
	for(uint32_t i=0;i<count;i++) {
		for(uint16_t j=0;j<dimension;j++)
			out << vectors[i](j) << ",";
		out << endl;
	}
}

