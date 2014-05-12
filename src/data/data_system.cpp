#include <DRDSP/data/data_system.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace DRDSP;

DataSystem::DataSystem() :
	dimension( 0 ),
	numParameters( 0 ),
	parameterDimension( 0 )
{}

DataSystem::DataSystem( uint32_t dim, uint32_t numParams, uint32_t paramDim ) :
	dimension( dim ),
	numParameters( numParams ),
	parameterDimension( paramDim ),
	dataSets( numParams ),
	parameters( numParams )
{}

//! Load a data system from a text file
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
	
	parameterDimension = (uint32_t)pDim;
	cout << "State Dimension: " << dimension << endl;
	cout << "Parameter Samples: " << numParameters << endl;
	cout << "Parameter Dimension: " << pDim << endl << endl;

	dataSets.resize( numParameters );
	parameters.resize( numParameters );

	string setPath;
	for(uint16_t i=0;i<numParameters;++i) {
		in >> setPath;
		if(binary) LoadSetBinary(setPath.c_str(),i,maxPoints);
		else LoadSetText(setPath.c_str(),i,maxPoints);
	}

	return true;
}

//! Load the ith data set from the given file in binary format
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
		++k;
	}
	return true;
}

//! Load the ith data set from the given file in text format (space separated)
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

	parameters[i].setZero( parameterDimension );
	for(uint16_t j=0;j<parameterDimension;++j) {
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

//! Write the data sets to numbered files with given prefix and suffix
void DataSystem::WriteDataSetsCSV( const char* filePrefix, const char* fileSuffix ) const {
	stringstream name;
	for(uint16_t i=0;i<numParameters;++i) {
		name.str("");
		name << filePrefix << parameters[i](0) << fileSuffix;
		dataSets[i].WriteCSV(name.str().c_str());
	}
}

//! Apply the given projection to this data system
DataSystem DataSystem::ProjectData( const MatrixXd& W ) const {
	DataSystem projectedData( (uint32_t)W.cols(), numParameters, parameterDimension );

	projectedData.parameters = parameters;

	for(uint16_t i=0;i<numParameters;++i) {
		projectedData.dataSets[i] = dataSets[i].ProjectData( W );
	}

	return projectedData;
}