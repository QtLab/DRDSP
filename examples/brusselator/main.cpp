#include <iostream>
#include <DRDSP/data/data_set.h>
#include <DRDSP/data/secants.h>
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/model_reduced_producer.h>

#include "brusselator.h"

using namespace std;
using namespace DRDSP;

struct Options {
	Options() : numIterations(10), maxPoints(20), targetDimension(2), numRBFs(30) {}
	uint32_t numIterations, maxPoints;
	uint16_t targetDimension, numRBFs;
};

Options GetOptions( int argc, char** argv ) {
	Options options;

	if( argc >= 2 ) options.targetDimension = (uint16_t)atoi(argv[1]);
	if( argc >= 3 ) options.numRBFs = (uint16_t)atoi(argv[2]);
	if( argc >= 4 ) options.maxPoints = (uint32_t)atoi(argv[3]);
	if( argc >= 5 ) options.numIterations = (uint32_t)atoi(argv[4]);

	return options;
}

void Compare( const ReducedDataSystem& reducedData, const DataSystem& rdata ) {

	ofstream out("output/comparison.csv");
	out << "Parameter,RMS,Max,Differences" << endl;
	for(uint16_t i=0;i<reducedData.numParameters;++i) {
		DataComparisonResult r = CompareData( reducedData.reducedData[i].points, rdata.dataSets[i].points );
		cout << "Parameter " << rdata.parameters[i] << endl;
		cout << "RMS: " << r.rmsDifference << endl;
		cout << "Max: " << r.maxDifference << endl;
		
		out << rdata.parameters[i] << ",";
		out << r.rmsDifference << ",";
		out << r.maxDifference << ",";
		for( const auto& x : r.differences )
			out << x << ",";
		out << endl;
	}

}

int main( int argc, char** argv ) {

	Options options = GetOptions(argc,argv);

	// Create a new data set
	DataSystem data;
	data.maxPoints = options.maxPoints;

	// Load data from file
	bool success = data.Load("data/brusselator.txt");
	if( !success ) {
		return 0;
	}

	// Pre-compute secants
	SecantsPreComputed* secants = new SecantsPreComputed [data.numParameters];
	SecantsPreComputed* newSecants = new SecantsPreComputed [data.numParameters];
	for(uint16_t i=0;i<data.numParameters;i++)
		secants[i].ComputeFromData( data.dataSets[i] );

	// Secant culling
	for(uint16_t i=0;i<data.numParameters;i++)
		newSecants[i] = secants[i].CullSecantsDegrees( 10.0 );

	delete[] secants;

	// Find a projection
	ProjSecant projSecant;
	projSecant.targetDimension = options.targetDimension;
	projSecant.targetMinProjectedLength = 0.7;

	// Compute initial condition
	projSecant.GetInitial(data);

	// Optimize over Grassmannian
	projSecant.Find(newSecants,data.numParameters);

	// Print some statistics
	projSecant.AnalyseSecants(newSecants,data.numParameters);

	delete[] newSecants;

	// Dynamics
	Brusselator brusselator;

	// Compute projected data
	cout << endl << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData(brusselator,data,projSecant.W);

	// Obtain the reduced model
	cout << endl << "Computing Reduced Model..." << endl;
	ModelReducedProducer producer;
	producer.numRBFs = options.numRBFs;
	ModelReduced reducedModel = producer.BruteForce(reducedData,data.parameterDimension,data.parameters.data(),options.numIterations);

	cout << "Total Cost = " << producer.ComputeTotalCost(reducedModel,reducedData,data.parameters.data()) << endl;

	reducedModel.WriteCSV("output/reduced.csv");
	//reducedData.WritePointsText("output/p2.5-points.csv");
	//reducedData.WriteVectorsText("output/p2.5-vectors.csv");
	projSecant.WriteBinary("output/projection.bin");
	projSecant.WriteCSV("output/projection.csv");

	//Compare( reducedData, rdata );

	return 0;
}
