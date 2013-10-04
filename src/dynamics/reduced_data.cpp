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

ReducedData::ReducedData() : points(nullptr), vectors(nullptr), derivatives(nullptr), dimension(0), count(0) {
}

ReducedData::ReducedData( uint32_t dim, uint32_t numPoints ) : points(nullptr), vectors(nullptr), derivatives(nullptr), dimension(0), count(0) {
	Create(dim,numPoints);
}

ReducedData::ReducedData( const ReducedData& rhs ) : points(nullptr), vectors(nullptr), derivatives(nullptr), dimension(0), count(0) {
	*this = rhs;			
}

ReducedData::ReducedData( ReducedData&& rhs ) {
	points = rhs.points;
	vectors = rhs.vectors;
	derivatives = rhs.derivatives;
	count = rhs.count;
	dimension = rhs.dimension;
	rhs.points = nullptr;
	rhs.vectors = nullptr;
	rhs.derivatives = nullptr;
	rhs.count = 0;
	rhs.dimension = 0;		
}

ReducedData::~ReducedData() {
	Destroy();
}

ReducedData& ReducedData::operator=( const ReducedData& rhs ) {
	Create(rhs.dimension,rhs.count);
	for(uint32_t i=0;i<count;i++) {
		points[i] = rhs.points[i];
		vectors[i] = rhs.vectors[i];
		derivatives[i] = rhs.derivatives[i];
	}
	return *this;
}

ReducedData& ReducedData::operator=( ReducedData&& rhs ) {
	if( this != &rhs ) {
		Destroy();
		points = rhs.points;
		vectors = rhs.vectors;
		derivatives = rhs.derivatives;
		count = rhs.count;
		dimension = rhs.dimension;
		rhs.points = nullptr;
		rhs.vectors = nullptr;
		rhs.derivatives = nullptr;
		rhs.count = 0;
		rhs.dimension = 0;
	}
	return *this;
}

void ReducedData::Create( uint32_t dim, uint32_t numPoints ) {
	if( count != numPoints ) {
		Destroy();
		points = new VectorXd [numPoints];
		vectors = new VectorXd [numPoints];
		derivatives = new MatrixXd [numPoints];
		count = numPoints;
	}
	dimension = dim;
	for(uint32_t i=0;i<count;i++) {
		points[i].setZero(dimension);
		vectors[i].setZero(dimension);
		derivatives[i].setZero(dimension,dimension);
	}
}

void ReducedData::Destroy() {
	delete[] points;
	points = nullptr;
	delete[] vectors;
	vectors = nullptr;
	delete[] derivatives;
	derivatives = nullptr;
	count = 0;
}

void ReducedData::ComputeData( ModelOriginal& original, const DataSet& data, const VectorXd& parameter, const MatrixXd& W ) {
	Create( W.cols(), data.count );

	static double stabilityFactor = 1.0;

	MatrixXd Ar, ArArT, I;

	for(uint32_t i=0;i<count;i++) {
		original.PrepareOptimizations( data.points[i], parameter );

		points[i] = W.adjoint() * data.points[i];
		vectors[i] = W.adjoint() * original.VectorField( data.points[i], parameter );

		if( original.componentWise ) {
			if( original.Euclidean ) {
				for(uint32_t c1=0;c1<original.oDim;c1++)
					for(uint32_t c2=0;c2<original.oDim;c2++) {
						double temp = original.VectorFieldDCW( data.points[i], parameter, c1, c2 );
						if( temp != 0.0 ) {
							for(uint32_t k=0;k<dimension;k++)
								for(uint32_t j=0;j<dimension;j++)
									derivatives[i](j,k) += W.adjoint()(j,c1) * temp * W(c2,k);
						}
					}
			}
		} else {
			if( original.Euclidean ) {
				derivatives[i] = W.adjoint() * original.VectorFieldD( data.points[i], parameter ) * W;
			} else {
				JacobiSVD<MatrixXd> svd(W.adjoint() * original.embedding->Derivative(original.embedding->Inverse(data.points[i])),ComputeThinU);
				uint32_t rank = 0;
				double tolerance = original.oDim * eps(svd.singularValues()(0));
				for(int j=0;j<svd.nonzeroSingularValues();j++) {
					if( svd.singularValues()(j) > tolerance )
						rank++;
					else break;
				}
				Ar = svd.matrixU().leftCols(rank) * svd.matrixU().leftCols(rank).transpose();
				I.setIdentity(Ar.rows(),Ar.cols());
				derivatives[i] = W.adjoint() * original.VectorFieldD( data.points[i], parameter ) * W * Ar - stabilityFactor*(I-Ar);
			}

		}

	}
}

void ReducedData::ComputeBoundingBox( VectorXd& bMin, VectorXd& bMax ) const {
	static double bboxScale = 1.0;
	bMin.setZero(dimension);
	bMax.setZero(dimension);
	for(uint32_t k=0;k<count;k++) {
		for(uint32_t j=0;j<dimension;j++) {
			if( points[k](j) > bMax(j) )
				bMax(j) = points[k](j);
			if( points[k](j) < bMin(j) )
				bMin(j) = points[k](j);
		}
	}
	// Scale it a bit
	VectorXd middle = (bMax + bMin) * 0.5;
	VectorXd halfDiff = (bMax - bMin) * 0.5;
	bMin = middle - halfDiff * bboxScale;
	bMax = middle + halfDiff * bboxScale;
}

bool ReducedData::ReadData( const char* filename ) {
	ifstream in(filename,ios::binary);
	if( !in ) {
		cout << "ReducedData::ReadData : file error" << endl;
		return false;
	}

	in.seekg(0, ios::end);
	if( in.tellg() < sizeof(double) * count * dimension * ( dimension + 2 ) ) {
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
	in.close();
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
	out.close();
}

