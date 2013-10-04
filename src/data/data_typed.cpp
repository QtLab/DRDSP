#include <DRDSP/data/data_typed.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace DRDSP;

DataSetTyped::DataSetTyped() {
}

void DataSetTyped::Create(uint32_t numPoints, uint32_t dim) {
	DataSet::Create(numPoints,dim);
	types = new uint16_t [count];
	for(uint32_t i=0;i<count;i++)
		types[i] = 0;
}

void DataSetTyped::Destroy() {
	DataSet::Destroy();
	delete[] types;
}

DataSetTyped::~DataSetTyped() {
	Destroy();
}

bool DataSystemTyped::Load( const char* filename ) {
	ifstream in(filename);
	if( !in )
		return false;
	numParameters = 0;
	uint16_t dt;
	bool binary;
	in >> binary >> dt >> numParameters >> dimension >> parameterDimension;
	if( !numParameters )
		return false;
	cout << "Loading system " << filename << endl << endl;
	dataSets = new DataSetTyped[numParameters];
	parameters = new VectorXd[numParameters];

	char* setPath = new char[256];
	stringstream fn;
	for(uint32_t i=0;i<numParameters;i++) {
		in >> setPath;
		fn.str("");
		fn << "data/" << setPath;
		ifstream inset(setPath);

		if(binary) LoadSetBinary(setPath,i);
		else LoadSetText(setPath,i);
	}
	delete[] setPath;
	in.close();
	return true;
}

bool DataSystemTyped::LoadSetBinary( const char* filename, uint32_t i ) {
	ifstream in;
	in.open(filename,ios::binary);
	if( !in )
		return false;

	cout << "Loading data file " << filename << endl << endl;

	uint32_t numPoints;
	in.read((char*)&numPoints,sizeof(numPoints));

	if( maxPoints && numPoints > maxPoints ) numPoints = maxPoints;
	
	parameters[i].setZero(parameterDimension);
	for( uint32_t j=0;j<parameterDimension;j++) {
		in.read((char*)&parameters[i](j),sizeof(double));
	}
	
	dataSets[i].Create(numPoints,dimension);
	uint32_t k=0, j=0;
	while( !in.eof() ) {
		if( k >= dataSets[i].count ) break;
		j=0;
		while( !in.eof() ) {
			if( j >= dimension ) break;
			in.read((char*)&dataSets[i].points[k](j++),sizeof(double));
		}
		in.read((char*)&dataSets[i].types[k],sizeof(uint8_t));
		k++;
	}
	in.close();
	return true;
}

bool DataSystemTyped::LoadSetText( const char* filename, uint32_t i ) {
	ifstream in;
	in.open(filename);
	if( !in )
		return false;

	cout << "Loading data file " << filename << endl << endl;

	uint32_t numPoints;
	in >> numPoints;

	if( maxPoints && numPoints > maxPoints ) numPoints = maxPoints;
	
	parameters[i].setZero(parameterDimension);
	for( uint32_t j=0;j<parameterDimension;j++) {
		in >> parameters[i](j);
	}
	
	dataSets[i].Create(numPoints,dimension);
	uint32_t k=0, j=0;
	while( !in.eof() ) {
		if ( k >= dataSets[i].count ) break;
		j=0;
		while( !in.eof() ) {
			if( j >= dimension ) break;
			in >> dataSets[i].points[k](j++);
		}
		in >> dataSets[i].types[k];
		k++;
	}
	in.close();
	return true;
}



