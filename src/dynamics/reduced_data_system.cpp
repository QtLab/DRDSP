#include <DRDSP/dynamics/reduced_data_system.h>

using namespace DRDSP;

ReducedDataSystem::ReducedDataSystem() : reducedData(nullptr), numParameters(0) {
}

ReducedDataSystem::ReducedDataSystem( uint32_t N ) : reducedData(nullptr), numParameters(0) {
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
	for(uint32_t i=0;i<numParameters;i++) {
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

void ReducedDataSystem::Create( uint32_t N ) {
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

void ReducedDataSystem::ComputeData( ModelOriginal& original, const DataSystem& data, const MatrixXd& W ) {
	Create( data.numParameters );
	for(uint32_t i=0;i<numParameters;i++) {
		reducedData[i].ComputeData(original,data.dataSets[i],data.parameters[i],W);
	}
}

void ReducedDataSystem::ComputeBoundingBox( VectorXd& bMin, VectorXd& bMax ) const {
	reducedData[0].ComputeBoundingBox(bMin,bMax);
	VectorXd bMin2, bMax2;
	uint32_t d = bMin.size();
	for(uint32_t i=1;i<numParameters;i++) {
		reducedData[i].ComputeBoundingBox(bMin2,bMax2);
		for(uint32_t j=0;j<d;j++) {
			if( bMin2(j) < bMin(j) ) bMin(j) = bMin2(j);
			if( bMax2(j) > bMax(j) ) bMax(j) = bMax2(j);
		}
	}
}

