#include <DRDSP/dynamics/reduced_data_system.h>
#include <DRDSP/dynamics/monodromy.h>
#include <sstream>
#include <fstream>
#include <iostream>

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

const ReducedDataSystem& ReducedDataSystem::WritePointsCSV( const char* filePrefix, const char* fileSuffix ) const {
	stringstream name;
	for(uint32_t i=0;i<numParameters;++i) {
		name.str("");
		name << filePrefix << i << fileSuffix;
		reducedData[i].WritePointsCSV(name.str().c_str());
	}
	return *this;
}

const ReducedDataSystem& ReducedDataSystem::WriteVectorsCSV( const char* filePrefix, const char* fileSuffix ) const {
	stringstream name;
	for(uint32_t i=0;i<numParameters;++i) {
		name.str("");
		name << filePrefix << i << fileSuffix;
		reducedData[i].WriteVectorsCSV(name.str().c_str());
	}
	return *this;
}

const ReducedDataSystem& ReducedDataSystem::WriteDerivativesCSV( const char* filePrefix, const char* fileSuffix ) const {
	stringstream name;
	for(uint32_t i=0;i<numParameters;++i) {
		name.str("");
		name << filePrefix << i << fileSuffix;
		reducedData[i].WriteDerivativesCSV(name.str().c_str());
	}
	return *this;
}

size_t ReducedDataSystem::TotalPoints() const {
	size_t N = 0;
	for(uint32_t i=0;i<numParameters;++i) {
		N += reducedData[i].count;
	}
	return N;
}

void DRDSP::Compare( const ReducedDataSystem& reducedData, const DataSystem& rdata ) {

	ofstream out("output/comparison.csv");
	out << "Parameter,RMS,Max,MaxMin,Scale,RMS,Max,MaxMin,Differences" << endl;
	for(uint32_t i=0;i<reducedData.numParameters;++i) {
		BallXd ball = ComputeBall( reducedData.reducedData[i].points );
		DataComparisonResult r = CompareData( reducedData.reducedData[i].points, rdata.dataSets[i].points );
		cout << "Parameter " << rdata.parameters[i] << endl;
		cout << "Max Rel Error: " << r.maxDifference / ball.radius << endl;

		out << rdata.parameters[i] << ",";
		out << r.rmsDifference << ",";
		out << r.maxDifference << ",";
		out << r.maxMinDifference << ",";
		out << ball.radius << ",";
		out << r.rmsDifference / ball.radius << ",";
		out << r.maxDifference / ball.radius << ",";
		out << r.maxMinDifference / ball.radius << ",";
		for( const auto& x : r.differences )
			out << x << ",";
		out << ",";
		for( const auto& x : r.differences )
			out << x / ball.radius << ",";
		out << endl;
	}
}

void DRDSP::ComparePeriods( const ReducedDataSystem& reducedData, const DataSystem& rdata, double dt, double tolerance ) {
	ofstream out("output/periods.csv");
	out << "Parameter,Original Period,Reduced Period,Error,Rel Error,%" << endl;
	for(uint32_t i=0;i<reducedData.numParameters;++i) {
		double p1 = DetectPeriod( reducedData[i].points, reducedData[i].vectors[0], dt, tolerance );
		double p2 = DetectPeriod( rdata.dataSets[i].points, reducedData[i].vectors[0], dt, tolerance );
		double err = p2 - p1;
		double rel = err / p1;
		double pc = rel * 100.0;
		out << rdata.parameters[i] << ",";
		out << p1 << ",";
		out << p2 << ",";
		out << err << ",";
		out << rel << ",";
		out << pc << ",";
		out << endl;
	}
}
