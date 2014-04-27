#include <DRDSP/dynamics/reduced_data_system.h>
#include <sstream>

using namespace DRDSP;
using namespace std;

ReducedDataSystem::ReducedDataSystem() : numParameters(0) {}

ReducedDataSystem::ReducedDataSystem( uint32_t N ) : reducedData(N), numParameters(N) {}

void ReducedDataSystem::Create( uint32_t N ) {
	reducedData.resize(N);
	numParameters = N;
}

AABB ReducedDataSystem::ComputeBoundingBox() const {
	AABB box = reducedData[0].ComputeBoundingBox();
	uint32_t d = reducedData[0].dimension;
	for(uint32_t i=1;i<numParameters;++i) {
		AABB box2 = reducedData[i].ComputeBoundingBox();
		for(uint32_t j=0;j<d;++j) {
			if( box2.bMin(j) < box.bMin(j) ) box.bMin(j) = box2.bMin(j);
			if( box2.bMax(j) > box.bMax(j) ) box.bMax(j) = box2.bMax(j);
		}
	}
	return box;
}

void ReducedDataSystem::WritePointsCSV( const char* filePrefix, const char* fileSuffix ) const {
	stringstream name;
	for(uint32_t i=0;i<numParameters;++i) {
		name.str("");
		name << filePrefix << i << fileSuffix;
		reducedData[i].WritePointsCSV(name.str().c_str());
	}
}

void ReducedDataSystem::WriteVectorsCSV( const char* filePrefix, const char* fileSuffix ) const {
	stringstream name;
	for(uint32_t i=0;i<numParameters;++i) {
		name.str("");
		name << filePrefix << i << fileSuffix;
		reducedData[i].WriteVectorsCSV(name.str().c_str());
	}
}

