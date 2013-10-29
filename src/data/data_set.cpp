#include <DRDSP/data/data_set.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace DRDSP;

DataSet::DataSet() : count(0), dimension(0), points(nullptr) {
}

DataSet::DataSet( const DataSet& rhs ) : count(0), dimension(0), points(nullptr) {
	Create(rhs.count,rhs.dimension);
	for(uint32_t i=0;i<count;i++)
		points[i] = rhs.points[i];
}

DataSet::DataSet( DataSet&& rhs ) {
	points = rhs.points;
	count = rhs.count;
	dimension = rhs.dimension;
	rhs.points = nullptr;
	rhs.count = 0;
}

DataSet::~DataSet() {
	Destroy();
}

DataSet& DataSet::operator=( const DataSet& rhs ) {
	Destroy();
	Create(rhs.count,rhs.dimension);
	for(uint32_t i=0;i<count;i++)
		points[i] = rhs.points[i];
	return *this;
}

DataSet& DataSet::operator=( DataSet&& rhs ) {
	if( this != &rhs ) {
		Destroy();
		points = rhs.points;
		count = rhs.count;
		dimension = rhs.dimension;
		rhs.points = nullptr;
		rhs.count = 0;
	}
	return *this;
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

	projectedData.Create(count,(uint32_t)W.cols());

	for(uint32_t i=0;i<count;i++) {
		projectedData.points[i] = W.adjoint() * points[i];
	}

	return std::move(projectedData);
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

	Create(numPoints,dimension);
	uint32_t k=0, j=0;
	while( !in.eof() ) {
		if( k >= count ) break;
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

	Create(numPoints,dimension);
	uint32_t k=0, j=0;
	while( !in.eof() ) {
		if ( k >= count ) break;
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

void DataSet::WriteText( const char* filename ) const {
	ofstream out(filename);
	if( !out ) {
		cout << "DataSet::WriteText : file error" << endl;
		return;
	}
	for(uint32_t i=0;i<count;i++) {
		for(uint16_t j=0;j<dimension;j++)
			out << points[i](j) << ",";
		out << endl;
	}
	out.close();
}

DataSystem::DataSystem() : dimension(0), numParameters(0), parameterDimension(0), maxPoints(0), dataSets(nullptr), parameters(nullptr) {
}

DataSystem::DataSystem( const DataSystem& rhs ) {
	Create(rhs.dimension,rhs.numParameters,rhs.parameterDimension);
	maxPoints = rhs.maxPoints;
	for(uint32_t i=0;i<numParameters;i++) {
		dataSets[i] = rhs.dataSets[i];
		parameters[i] = rhs.parameters[i];
	}
}

DataSystem::DataSystem( DataSystem&& rhs ) {
	dataSets = rhs.dataSets;
	parameters = rhs.parameters;
	dimension = rhs.dimension;
	numParameters = rhs.numParameters;
	parameterDimension = rhs.parameterDimension;
	maxPoints = rhs.maxPoints;
	rhs.dataSets = nullptr;
	rhs.parameters = nullptr;
	rhs.numParameters = 0;
}

DataSystem::~DataSystem() {
	Destroy();
}

DataSystem& DataSystem::operator=( const DataSystem& rhs ) {
	Destroy();
	Create(rhs.dimension,rhs.numParameters,rhs.parameterDimension);
	maxPoints = rhs.maxPoints;
	for(uint32_t i=0;i<numParameters;i++) {
		dataSets[i] = rhs.dataSets[i];
		parameters[i] = rhs.parameters[i];
	}
	return *this;
}

DataSystem& DataSystem::operator=( DataSystem&& rhs ) {
	if( this != &rhs ) {
		Destroy();
		dataSets = rhs.dataSets;
		parameters = rhs.parameters;
		dimension = rhs.dimension;
		numParameters = rhs.numParameters;
		parameterDimension = rhs.parameterDimension;
		maxPoints = rhs.maxPoints;
		rhs.dataSets = nullptr;
		rhs.parameters = nullptr;
		rhs.numParameters = 0;
	}
	return *this;
}

void DataSystem::Create( uint32_t dim, uint16_t numParams, uint8_t paramDim ) {
	numParameters = numParams;
	dimension = dim;
	parameterDimension = paramDim;
	dataSets = new DataSet [numParameters];
	parameters = new VectorXd [numParameters];
	for(uint16_t i=0;i<numParameters;i++) {
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

	dataSets = new DataSet [numParameters];
	parameters = new VectorXd [numParameters];

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
	ifstream in;
	in.open(filename,ios::binary);
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
	
	dataSets[i].Create(numPoints,dimension);
	uint32_t k=0;
	while( !in.eof() ) {
		if( k >= dataSets[i].count ) break;
		in.read((char*)&dataSets[i].points[k](0),sizeof(double)*dimension);
		k++;
	}
	in.close();
	return true;
}

bool DataSystem::LoadSetText( const char* filename, uint32_t i ) {
	ifstream in;
	in.open(filename);
	if( !in ) {
		cout << "File not found" << endl;
		return false;
	}

	cout << "Loading data file " << filename;

	uint32_t numPoints;
	in >> numPoints;

	if( maxPoints && numPoints > maxPoints ) numPoints = maxPoints;

	parameters[i].setZero(parameterDimension);
	for(uint16_t j=0;j<parameterDimension;j++) {
		in >> parameters[i](j);
	}

	cout << " - points = " << numPoints << ", parameter(0) = " << parameters[i](0) << endl << endl;

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

	projectedData.Create((uint32_t)W.cols(),numParameters,parameterDimension);

	for(uint16_t i=0;i<numParameters;i++) {
		projectedData.parameters[i] = parameters[i];
		projectedData.dataSets[i].Create(dataSets[i].count,projectedData.dimension);

		for(uint32_t j=0;j<dataSets[i].count;j++) {
			projectedData.dataSets[i].points[j] = W.adjoint() * dataSets[i].points[j];
		}

	}

	return std::move(projectedData);
}