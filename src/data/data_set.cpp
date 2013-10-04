#include <DRDSP/data/data_set.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace DRDSP;

DataSet::DataSet() : count(0), dimension(0) {
}

DataSet::DataSet( const DataSet& rhs ) {
	Create(rhs.count,rhs.dimension);
	for(uint32_t i=0;i<count;i++)
		points[i] = rhs.points[i];
}

DataSet::DataSet( DataSet&& rhs ) {
	delete[] points;
	points = rhs.points;
	count = rhs.count;
	dimension = rhs.dimension;
	rhs.points = nullptr;
	rhs.count = 0;
}

DataSet::~DataSet() {
	Destroy();
}

void DataSet::Create(uint32_t numPoints, uint32_t dim) {
	count = numPoints;
	dimension = dim;
	points = new VectorXd [count];
	for(uint32_t i=0;i<count;i++)
		points[i].setZero(dimension);
}

void DataSet::Destroy() {
	delete[] points;
	points = nullptr;
	count = 0;
}

DataSet DataSet::ProjectData( const MatrixXd& W ) const {
	DataSet projectedData;

	projectedData.Create(count,W.cols());

	for(uint32_t i=0;i<count;i++) {
		projectedData.points[i] = W.adjoint() * points[i];
	}

	return std::move(projectedData);
}

DataSystem::DataSystem() : numParameters(0), dimension(0), parameterDimension(0) {
}

DataSystem::DataSystem( const DataSystem& rhs ) {
	Create(rhs.numParameters,rhs.parameterDimension,rhs.dimension);
	for(uint32_t i=0;i<numParameters;i++)
		dataSets[i] = rhs.dataSets[i];
}

DataSystem::DataSystem( DataSystem&& rhs ) {
	delete[] dataSets;
	delete[] parameters;
	dataSets = rhs.dataSets;
	parameters = rhs.parameters;
	dimension = rhs.dimension;
	numParameters = rhs.numParameters;
	parameterDimension = rhs.parameterDimension;
	rhs.dataSets = nullptr;
	rhs.parameters = nullptr;
	rhs.numParameters = 0;
}

DataSystem::~DataSystem() {
	Destroy();
}

void DataSystem::Create( uint32_t numParams, uint32_t paramDim, uint32_t dim ) {
	numParameters = numParams;
	dimension = dim;
	parameterDimension = paramDim;
	dataSets = new DataSet [numParameters];
	parameters = new VectorXd [numParameters];
	for(uint32_t i=0;i<numParameters;i++) {
		parameters[i].setZero(parameterDimension);
	}

}

void DataSystem::Destroy() {
	delete[] dataSets;
	delete[] parameters;
	dataSets = nullptr;
	parameters = nullptr;
	numParameters = 0;
}

bool DataSystem::Load( const char* filename ) {
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
	dataSets = new DataSet[numParameters];
	parameters = new VectorXd[numParameters];

	char* setPath = new char[256];
	stringstream fn;
	for(uint32_t i=0;i<numParameters;i++) {
		in >> setPath;
		fn.str("");
		fn << "data/" << setPath;

		if(binary) LoadSetBinary(setPath,i);
		else LoadSetText(setPath,i);
	}
	delete[] setPath;
	in.close();
	return true;
}


bool DataSystem::LoadSetBinary( const char* filename, uint32_t i ) {
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
		k++;
	}
	in.close();
	return true;
}

bool DataSystem::LoadSetText( const char* filename, uint32_t i ) {
	ifstream in;
	in.open(filename);
	if( !in )
		return false;

	cout << "Loading data file " << filename << endl << endl;

	uint32_t numPoints;
	in >> numPoints;

	if( maxPoints && numPoints > maxPoints ) numPoints = maxPoints;
	
	parameters[i].setZero(parameterDimension);
	for(uint32_t j=0;j<parameterDimension;j++) {
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
		k++;
	}
	in.close();
	return true;
}

DataSystem DataSystem::ProjectData( const MatrixXd& W ) const {
	DataSystem projectedData;

	projectedData.Create(numParameters,parameterDimension,W.cols());

	for(uint32_t i=0;i<numParameters;i++) {
		projectedData.parameters[i] = parameters[i];
		projectedData.dataSets[i].Create(dataSets[i].count,projectedData.dimension);

		for(uint32_t j=0;j<dataSets[i].count;j++) {
			projectedData.dataSets[i].points[j] = W.adjoint() * dataSets[i].points[j];
		}

	}

	return std::move(projectedData);
}


