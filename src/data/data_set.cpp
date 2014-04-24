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

DataSet DataSet::ProjectData( const MatrixXd& W ) const {
	DataSet projectedData( points.size(), (uint32_t)W.cols() );

	for(uint32_t i=0;i<points.size();++i) {
		projectedData.points[i] = W.adjoint() * points[i];
	}

	return projectedData;
}

bool DataSet::LoadBinary( const char* filename ) {
	ifstream in;
	in.open(filename,ios::binary);
	if( !in ) {
		cout << "File not found: " << filename << endl;
		return false;
	}

	cout << "Loading data file: " << filename;

	uint32_t numPoints;
	in.read((char*)&numPoints,sizeof(numPoints));
	cout << " - " << numPoints << " points" << endl << endl;

	in.seekg(sizeof(double),ios::cur);

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
		k++;
	}
	in.close();
	return true;
}

bool DataSet::LoadText( const char* filename ) {
	ifstream in;
	in.open(filename);
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
		k++;
	}
	in.close();
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
	out.close();
}

DataSystem::DataSystem() :
	dimension( 0 ),
	numParameters( 0 ),
	parameterDimension( 0 ),
	maxPoints( 0 )
{}

DataSystem::DataSystem( uint32_t dim, uint16_t numParams, uint8_t paramDim ) :
	dimension( dim ),
	numParameters( numParams ),
	parameterDimension( paramDim ),
	maxPoints( 0 ),
	dataSets( numParameters ),
	parameters( numParameters )
{
	for( auto& p : parameters ) {
		p.setZero(parameterDimension);
	}
}

bool DataSystem::Load( const char* filename ) {
	ifstream in(filename);
	if( !in ) {
		cout << "File not found" << endl;
		return false;
	}
	cout << "Loading system " << filename << endl << endl;
	numParameters = 0;
	uint16_t dt, pDim;
	bool binary;
	in >> binary >> dt >> numParameters >> dimension >> pDim;
	if( !numParameters )
		return false;
	
	parameterDimension = (uint8_t)pDim;
	cout << "State Dimension: " << dimension << endl;
	cout << "Parameter Samples: " << numParameters << endl;
	cout << "Parameter Dimension: " << pDim << endl << endl;

	dataSets.resize( numParameters );
	parameters.resize( numParameters );

	char* setPath = new char [256];
	for(uint16_t i=0;i<numParameters;i++) {
		in >> setPath;
		if(binary) LoadSetBinary(setPath,i);
		else LoadSetText(setPath,i);
	}
	delete[] setPath;
	in.close();
	return true;
}


bool DataSystem::LoadSetBinary( const char* filename, uint32_t i ) {
	ifstream in(filename,ios::binary);
	if( !in ) {
		cout << "File not found" << endl;
		return false;
	}

	cout << "Loading data file " << filename;

	uint32_t numPoints;
	in.read((char*)&numPoints,sizeof(numPoints));

	if( maxPoints && numPoints > maxPoints ) numPoints = maxPoints;

	parameters[i].setZero(parameterDimension);
	in.read((char*)&parameters[i](0),sizeof(double)*parameterDimension);

	cout << " - points = " << numPoints << ", parameter(0) = " << parameters[i](0) << endl << endl;
	
	dataSets[i].points.resize( numPoints );
	
	for( auto& x : dataSets[i].points )
		x.setZero( dimension );

	uint32_t k=0;
	while( !in.eof() ) {
		if( k >= dataSets[i].points.size() ) break;
		in.read((char*)&dataSets[i].points[k](0),sizeof(double)*dimension);
		k++;
	}
	return true;
}

bool DataSystem::LoadSetText( const char* filename, uint32_t i ) {
	ifstream in(filename);
	if( !in ) {
		cout << "File not found" << endl;
		return false;
	}

	cout << "Loading data file " << filename;

	uint32_t numPoints;
	in >> numPoints;

	if( maxPoints && numPoints > maxPoints ) numPoints = maxPoints;

	parameters[i].setZero( parameterDimension );
	for(uint16_t j=0;j<parameterDimension;j++) {
		in >> parameters[i](j);
	}

	cout << " - points = " << numPoints << ", parameter(0) = " << parameters[i](0) << endl << endl;

	dataSets[i].points.resize( numPoints );
	for( auto& x : dataSets[i].points )
		x.setZero( dimension );

	uint32_t k=0, j=0;
	while( !in.eof() ) {
		if ( k >= dataSets[i].points.size() ) break;
		j=0;
		while( !in.eof() ) {
			if( j >= dimension ) break;
			in >> dataSets[i].points[k](j++);
		}
		k++;
	}
	return true;
}

void DataSystem::WriteDataSetsCSV( const char* filePrefix, const char* fileSuffix ) const {
	stringstream name;
	for(uint16_t i=0;i<numParameters;i++) {
		name.str("");
		name << filePrefix << parameters[i](0) << fileSuffix;
		dataSets[i].WriteCSV(name.str().c_str());
	}
}

DataSystem DataSystem::ProjectData( const MatrixXd& W ) const {
	DataSystem projectedData( (uint32_t)W.cols(), numParameters, parameterDimension );

	for(uint16_t i=0;i<numParameters;i++) {
		projectedData.parameters[i] = parameters[i];
		projectedData.dataSets[i] = dataSets[i].ProjectData( W );
	}

	return projectedData;
}

vector<double> DRDSP::NormDifferences( const vector<VectorXd>& s1, const vector<VectorXd>& s2 ) {
	size_t n = min( s1.size(), s2.size() );
	vector<double> r(n);

	for(size_t i=0;i<n;++i) {
		r[i] = ( s1[i] - s2[i] ).norm();
	}

	return r;
}

double DRDSP::RMSDifference( const vector<VectorXd>& s1, const vector<VectorXd>& s2 ) {
	size_t n = min( s1.size(), s2.size() );
	double sum = 0.0;

	for(size_t i=0;i<n;++i) {
		sum += ( s1[i] - s2[i] ).squaredNorm();
	}

	return sqrt( sum / n );
}