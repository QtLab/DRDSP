#include <DRDSP/dynamics/reduced_data_system.h>

using namespace DRDSP;

ReducedDataSystem::ReducedDataSystem() : reducedData(nullptr), numParameters(0) {
}

ReducedDataSystem::ReducedDataSystem( uint16_t N ) : reducedData(nullptr), numParameters(0) {
	Create(N);
}

ReducedDataSystem::ReducedDataSystem( const ReducedDataSystem& rhs ) : reducedData(nullptr), numParameters(0) {
	*this = rhs;			
}

ReducedDataSystem::ReducedDataSystem( ReducedDataSystem&& rhs ) {
	reducedData = rhs.reducedData;
	numParameters = rhs.numParameters;
	rhs.reducedData = nullptr;
	rhs.numParameters = 0;		
}

ReducedDataSystem::~ReducedDataSystem() {
	Destroy();
}

ReducedDataSystem& ReducedDataSystem::operator=( const ReducedDataSystem& rhs ) {
	Create(rhs.numParameters);
	for(uint16_t i=0;i<numParameters;i++) {
		reducedData[i] = rhs.reducedData[i];
	}
	return *this;
}

ReducedDataSystem& ReducedDataSystem::operator=( ReducedDataSystem&& rhs ) {
	if( this != &rhs ) {
		Destroy();
		reducedData = rhs.reducedData;
		numParameters = rhs.numParameters;
		rhs.reducedData = nullptr;
		rhs.numParameters = 0;
	}
	return *this;
}

void ReducedDataSystem::Create( uint16_t N ) {
	if( numParameters != N ) {
		Destroy();
		reducedData = new ReducedData [N];
		numParameters = N;
	}
}

void ReducedDataSystem::Destroy() {
	delete[] reducedData;
	reducedData = nullptr;
	numParameters = 0;
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

