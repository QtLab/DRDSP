#include <fstream>
#include <iostream>
#include <DRDSP/dynamics/reduced_data.h>

using namespace std;
using namespace DRDSP;

double DRDSP::eps( double x ) {
	int64_t y = *((int64_t*)&x);
	y++;
	return *((double*)&y) - x;
}

ReducedData::ReducedData() : dimension(0), count(0) {}

ReducedData::ReducedData( uint32_t dim, size_t numPoints ) : dimension(0), count(0) {
	Create(dim,numPoints);
}

void ReducedData::Create( uint32_t dim, size_t numPoints ) {
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

AABB ReducedData::ComputeBoundingBox() const {
	AABB box(dimension);
	box.SetZero();
	for(uint32_t k=0;k<count;k++) {
		for(uint32_t j=0;j<dimension;j++) {
			if( points[k](j) > box.bMax(j) )
				box.bMax(j) = points[k](j);
			if( points[k](j) < box.bMin(j) )
				box.bMin(j) = points[k](j);
		}
	}
	return box;
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
	if( (size_t)in.tellg() < sizeof(double) * count * dimension * ( dimension + 2 ) ) {
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
		cout << "ReducedData::WritePointsCSV : file error" << endl;
		return;
	}
	for(uint32_t i=0;i<count;i++) {
		for(uint32_t j=0;j<dimension;j++)
			out << points[i](j) << ",";
		out << endl;
	}
}

void ReducedData::WriteVectorsCSV( const char* filename ) const {
	ofstream out(filename);
	if( !out ) {
		cout << "ReducedData::WriteVectorsCSV : file error" << endl;
		return;
	}
	for(uint32_t i=0;i<count;++i) {
		for(uint32_t j=0;j<dimension;++j)
			out << vectors[i](j) << ",";
		out << endl;
	}
}

void ReducedData::WriteDerivativesCSV( const char* filename ) const {
	ofstream out(filename);
	if( !out ) {
		cout << "ReducedData::WriteDerivativesCSV : file error" << endl;
		return;
	}
	for(uint32_t i=0;i<count;++i) {
		for(uint32_t j=0;j<dimension;++j)
			for(uint32_t k=0;k<dimension;++k)
				out << derivatives[i](j,k) << ",";
		out << endl;
	}
}
