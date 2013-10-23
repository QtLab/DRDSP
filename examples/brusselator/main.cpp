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

	// Pre-compute secants if less than 64 MB
	Secants* secants = new Secants [data.numParameters];
	Secants* newSecants = new Secants [data.numParameters];
	for(uint16_t i=0;i<data.numParameters;i++)
		secants[i].ComputeFromData( data.dataSets[i], 1 << 26 );

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
	ModelReduced reducedModel = producer.BruteForce(reducedData,data.parameterDimension,data.parameters,options.numIterations);

	cout << "Total Cost = " << producer.ComputeTotalCost(reducedModel,reducedData,data.parameters) << endl;

	reducedModel.OutputText("output/reduced.csv");
	//reducedData.WritePointsText("output/p2.5-points.csv");
	//reducedData.WriteVectorsText("output/p2.5-vectors.csv");
	projSecant.WriteBinary("output/projection.bin");
	projSecant.WriteText("output/projection.csv");

	return 0;
}
