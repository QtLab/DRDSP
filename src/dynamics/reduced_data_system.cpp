#include <DRDSP/dynamics/reduced_data_system.h>
#include <sstream>

using namespace DRDSP;
using namespace std;

ReducedDataSystem::ReducedDataSystem() : numParameters(0) {}

ReducedDataSystem::ReducedDataSystem( uint16_t N ) : reducedData(N), numParameters(N) {}

void ReducedDataSystem::Create( uint16_t N ) {
	reducedData.resize(N);
	numParameters = N;
}

void ReducedDataSystem::ComputeData( ModelParameterized& original, const DataSystem& data, const MatrixXd& W ) {
	Create( data.numParameters );
	for(uint16_t i=0;i<numParameters;i++) {
		reducedData[i].ComputeData(original,data.parameters[i],data.dataSets[i],W);
	}
}

void ReducedDataSystem::ComputeData( ModelParameterizedCW& original, const DataSystem& data, const MatrixXd& W ) {
	Create( data.numParameters );
	for(uint16_t i=0;i<numParameters;i++) {
		reducedData[i].ComputeData(original,data.parameters[i],data.dataSets[i],W);
	}
}

void ReducedDataSystem::ComputeData( ModelParameterizedEmbedded& original, const DataSystem& data, const MatrixXd& W ) {
	Create( data.numParameters );
	for(uint16_t i=0;i<numParameters;i++) {
		reducedData[i].ComputeData(original,data.parameters[i],data.dataSets[i],W);
	}
}

void ReducedDataSystem::ComputeData( ModelParameterizedEmbeddedCW& original, const DataSystem& data, const MatrixXd& W ) {
	Create( data.numParameters );
	for(uint16_t i=0;i<numParameters;i++) {
		reducedData[i].ComputeData(original,data.parameters[i],data.dataSets[i],W);
	}
}

AABB ReducedDataSystem::ComputeBoundingBox() const {
	AABB box = reducedData[0].ComputeBoundingBox();
	uint32_t d = reducedData[0].dimension;
	for(uint16_t i=1;i<numParameters;i++) {
		AABB box2 = reducedData[i].ComputeBoundingBox();
		for(uint16_t j=0;j<d;j++) {
			if( box2.bMin(j) < box.bMin(j) ) box.bMin(j) = box2.bMin(j);
			if( box2.bMax(j) > box.bMax(j) ) box.bMax(j) = box2.bMax(j);
		}
	}
	return std::move(box);
}

void ReducedDataSystem::WritePointsCSV( const char* filePrefix, const char* fileSuffix ) const {
	stringstream name;
	for(uint16_t i=0;i<numParameters;i++) {
		name.str("");
		name << filePrefix << i << fileSuffix;
		reducedData[i].WritePointsCSV(name.str().c_str());
	}
}

void ReducedDataSystem::WriteVectorsCSV( const char* filePrefix, const char* fileSuffix ) const {
	stringstream name;
	for(uint16_t i=0;i<numParameters;i++) {
		name.str("");
		name << filePrefix << i << fileSuffix;
		reducedData[i].WriteVectorsCSV(name.str().c_str());
	}
}

