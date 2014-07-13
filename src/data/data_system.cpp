#include <DRDSP/data/data_system.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace DRDSP;

DataSystem::DataSystem() :
	dimension( 0 ),
	numParameters( 0 ),
	paramDim( 0 )
{}

DataSystem::DataSystem( uint32_t dim, uint32_t numParams, uint32_t paramDim ) :
	dimension( dim ),
	numParameters( numParams ),
	paramDim( paramDim ),
	dataSets( numParams ),
	parameters( numParams )
{}

bool DataSystem::Load( const char* filename, uint32_t maxPoints ) {
	ifstream in(filename);
	if( !in ) {
		cout << "DataSystem::Load - File not found: " << filename << endl;
		return false;
	}
	cout << "Loading system " << filename << endl << endl;
	numParameters = 0;
	uint16_t dt, pDim;
	bool binary;
	in >> binary >> dt >> numParameters >> dimension >> pDim;
	if( !numParameters )
		return false;
	
	paramDim = (uint32_t)pDim;
	cout << "State Dimension: " << dimension << endl;
	cout << "Parameter Samples: " << numParameters << endl;
	cout << "Parameter Dimension: " << pDim << endl << endl;

	dataSets.resize( numParameters );
	parameters.resize( numParameters );

	string setPath;
	for(uint32_t i=0;i<numParameters;++i) {
		in >> setPath;
		if(binary) LoadSetBinary(setPath.c_str(),i,maxPoints);
		else LoadSetText(setPath.c_str(),i,maxPoints);
	}

	return true;
}

bool DataSystem::LoadSetBinary( const char* filename, uint32_t i, uint32_t maxPoints ) {
	ifstream in(filename,ios::binary);
	if( !in ) {
		cout << "File not found" << endl;
		return false;
	}

	cout << "Loading data file " << filename;

	uint32_t numPoints;
	in.read((char*)&numPoints,sizeof(numPoints));

	if( maxPoints && numPoints > maxPoints ) numPoints = maxPoints;

	parameters[i].setZero(paramDim);
	in.read((char*)&parameters[i](0),sizeof(double)*paramDim);

	cout << " - points = " << numPoints << ", parameter(0) = " << parameters[i](0) << endl << endl;
	
	dataSets[i].points.resize( numPoints );
	
	for( auto& x : dataSets[i].points )
		x.setZero( dimension );

	uint32_t k=0;
	while( !in.eof() ) {
		if( k >= dataSets[i].points.size() ) break;
		in.read((char*)&dataSets[i].points[k](0),sizeof(double)*dimension);
		++k;
	}
	return true;
}

bool DataSystem::LoadSetText( const char* filename, uint32_t i, uint32_t maxPoints ) {
	ifstream in(filename);
	if( !in ) {
		cout << "File not found" << endl;
		return false;
	}

	cout << "Loading data file " << filename;

	uint32_t numPoints;
	in >> numPoints;

	if( maxPoints && numPoints > maxPoints ) numPoints = maxPoints;

	parameters[i].setZero( paramDim );
	for(uint32_t j=0;j<paramDim;++j) {
		in >> parameters[i](j);
	}

	cout << " - points = " << numPoints << ", parameter(0) = " << parameters[i](0) << endl << endl;

	dataSets[i].points.resize( numPoints );
	for( auto& x : dataSets[i].points )
		x.setZero( dimension );

	uint32_t k = 0, j = 0;
	while( !in.eof() ) {
		if( k >= dataSets[i].points.size() ) break;
		j = 0;
		while( !in.eof() ) {
			if( j >= dimension ) break;
			in >> dataSets[i].points[k](j++);
		}
		++k;
	}
	return true;
}

void DataSystem::WriteCSV( const char* filePrefix, const char* fileSuffix ) const {
	stringstream name;
	for(uint32_t i=0;i<numParameters;++i) {
		name.str("");
		name << filePrefix << parameters[i](0) << fileSuffix;
		dataSets[i].WriteCSV(name.str().c_str());
	}
}

void DataSystem::WriteBinary( const char* filePrefix, const char* fileSuffix ) const {
	stringstream name;
	for(uint32_t i=0;i<numParameters;++i) {
		name.str("");
		name << filePrefix << parameters[i](0) << fileSuffix;
		dataSets[i].WriteBinary(name.str().c_str());
	}
}

DataSystem DataSystem::ProjectData( const MatrixXd& W ) const {
	DataSystem projectedData( (uint32_t)W.cols(), numParameters, paramDim );

	projectedData.parameters = parameters;

	for(uint32_t i=0;i<numParameters;++i) {
		projectedData.dataSets[i] = MapMatrix( W, dataSets[i] );
	}

	return projectedData;
}

size_t DataSystem::TotalPoints() const {
	size_t N = 0;
	for(uint32_t i=0;i<numParameters;++i) {
		N += dataSets[i].points.size();
	}
	return N;
}