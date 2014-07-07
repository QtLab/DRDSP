#include <DRDSP/data/data_set.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace DRDSP;

DataSet::DataSet( size_t numPoints, uint32_t dim ) : dimension(dim), points(numPoints) {
	for( auto& x : points )
		x.setZero(dimension);
}

bool DataSet::LoadBinary( const char* filename ) {
	ifstream in(filename,ios::binary);
	if( !in ) {
		cout << "File not found: " << filename << endl;
		return false;
	}

	cout << "Loading data file: " << filename;

	uint32_t numPoints;
	in.read((char*)&numPoints,sizeof(numPoints));
	cout << " - " << numPoints << " points" << endl << endl;

	//in.seekg(sizeof(double),ios::cur);

	points.resize( numPoints );
	for( auto& x : points )
		x.setZero( dimension );

	uint32_t k=0, j=0;
	while( !in.eof() ) {
		if( k >= points.size() ) break;
		j=0;
		while( !in.eof() ) {
			if( j >= dimension ) break;
			in.read((char*)&points[k](j++),sizeof(double));
		}
		++k;
	}
	return true;
}

bool DataSet::LoadText( const char* filename ) {
	ifstream in(filename);
	if( !in ) {
		cout << "File not found: " << filename << endl;
		return false;
	}

	cout << "Loading data file: " << filename;

	uint32_t numPoints;
	in >> numPoints;
	cout << " - " << numPoints << " points" << endl << endl;

	points.resize( numPoints );
	for( auto& x : points )
		x.setZero( dimension );

	uint32_t k=0, j=0;
	while( !in.eof() ) {
		if ( k >= points.size() ) break;
		j=0;
		while( !in.eof() ) {
			if( j >= dimension ) break;
			in >> points[k](j++);
		}
		++k;
	}
	return true;
}

void DataSet::WriteCSV( const char* filename ) const {
	ofstream out(filename);
	if( !out ) {
		cout << "DataSet::WriteText : file error" << endl;
		return;
	}
	for( const auto& x : points ) {
		for(uint32_t j=0;j<dimension;++j)
			out << x(j) << ",";
		out << endl;
	}
}

void DataSet::WriteBinary( const char* filename ) const {
	ofstream out(filename,ios::binary);
	if( !out ) {
		cout << "DataSet::WriteBinary : file error" << endl;
		return;
	}
	uint32_t numPoints = (uint32_t)points.size();
	out.write( (const char*)&numPoints, sizeof(uint32_t) );
	for( const auto& x : points ) {
		out.write( (const char*)&x[0], sizeof(double)*dimension );
	}
}

DataComparisonResult DRDSP::CompareData( const vector<VectorXd>& s1, const vector<VectorXd>& s2 ) {
	size_t n = min( s1.size(), s2.size() );
	DataComparisonResult r(n);
	double sum = 0.0, maxVal = 0.0;

	for(size_t i=0;i<n;++i) {
		double temp = ( s1[i] - s2[i] ).squaredNorm();
		r.differences[i] = sqrt( temp );
		sum += temp;
		if( r.differences[i] > maxVal ) maxVal = r.differences[i];
	}
	r.rmsDifference = sqrt( sum / n );
	r.maxDifference = maxVal;

	maxVal = 0.0;
	for( const auto& x : s1 ) {
		double minVal = numeric_limits<double>::max();
		for( const auto& y : s2 ) {
			double diff = ( x - y ).norm();
			if( diff < minVal ) minVal = diff;
		}
		if( minVal > maxVal ) maxVal = minVal;
	}
	r.maxMinDifference = maxVal;
	
	return r;
}

BallXd DRDSP::ComputeBall( const vector<VectorXd>& points ) {
	VectorXd sum(points[0].size());
	sum.setZero();
	for( const auto& x : points ) {
		sum += x;
	}
	BallXd ball;
	ball.centre = sum / double(points.size());
	double r2Max = 0.0;
	for( const auto& x : points ) {
		double r2 = ( x - ball.centre ).squaredNorm();
		if( r2 > r2Max ) r2Max = r2;
	}
	ball.radius = sqrt(r2Max);
	return ball;
}

