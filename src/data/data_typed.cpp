#include <DRDSP/data/data_typed.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace DRDSP;

DataSetTyped::DataSetTyped( uint32_t numPoints, uint32_t dim ) :
	DataSet( numPoints, dim ),
	types( numPoints, uint16_t(0) )
{}

DataSystemTyped::DataSystemTyped( uint32_t dim, uint16_t numParams, uint8_t paramDim ) :
	numParameters(numParams),
	dimension(dim),
	parameterDimension(paramDim),
	dataSets(numParams),
	parameters(numParams)
{
	for( auto& p : parameters ) {
		p.setZero(parameterDimension);
	}
}

bool DataSystemTyped::Load( const char* filename, uint32_t maxPoints ) {
	ifstream in(filename);
	if( !in ) {
		cout << "File not found" << endl;
		return false;
	}
	cout << "Loading system " << filename << endl << endl;
	numParameters = 0;
	uint16_t dt;
	bool binary;
	in >> binary >> dt >> numParameters >> dimension >> parameterDimension;
	if( !numParameters )
		return false;
	
	cout << "State Dimension: " << dimension << endl;
	cout << "Parameter Samples: " << numParameters << endl;
	cout << "Parameter Dimension: " << parameterDimension << endl << endl;

	dataSets.resize(numParameters);
	parameters.resize(numParameters);

	char* setPath = new char [256];
	for(uint16_t i=0;i<numParameters;i++) {
		in >> setPath;
		if(binary) LoadSetBinary(setPath,i);
		else LoadSetText(setPath,i);
	}
	delete[] setPath;
	return true;
}

bool DataSystemTyped::LoadSetBinary( const char* filename, uint32_t i, uint32_t maxPoints ) {
	ifstream in(filename,ios::binary);
	if( !in ) {
		cout << "File not found" << endl;
		return false;
	}

	cout << "Loading data file " << filename << endl << endl;

	uint32_t numPoints;
	in.read((char*)&numPoints,sizeof(numPoints));

	if( maxPoints && numPoints > maxPoints ) numPoints = maxPoints;
	
	parameters[i].setZero(parameterDimension);
	for(uint16_t j=0;j<parameterDimension;j++) {
		in.read((char*)&parameters[i](j),sizeof(double));
	}
	
	dataSets[i].points.resize( numPoints );
	for( auto& x : dataSets[i].points )
		x.setZero( dimension );

	uint32_t k=0, j=0;
	while( !in.eof() ) {
		if( k >= dataSets[i].points.size() ) break;
		j=0;
		while( !in.eof() ) {
			if( j >= dimension ) break;
			in.read((char*)&dataSets[i].points[k](j++),sizeof(double));
		}
		in.read((char*)&dataSets[i].types[k],sizeof(uint8_t));
		k++;
	}
	return true;
}

bool DataSystemTyped::LoadSetText( const char* filename, uint32_t i, uint32_t maxPoints ) {
	ifstream in(filename);
	if( !in ) {
		cout << "File not found" << endl;
		return false;
	}

	cout << "Loading data file " << filename << endl << endl;

	uint32_t numPoints;
	in >> numPoints;

	if( maxPoints && numPoints > maxPoints ) numPoints = maxPoints;
	
	parameters[i].setZero(parameterDimension);
	for(uint16_t j=0;j<parameterDimension;j++) {
		in >> parameters[i](j);
	}
	
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
		in >> dataSets[i].types[k];
		k++;
	}
	return true;
}



