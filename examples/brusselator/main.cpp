#include <iostream>
#include <DRDSP/data/data_set.h>
#include <DRDSP/data/secants.h>
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/model_reduced_producer.h>

#include "brusselator.h"

using namespace std;
using namespace DRDSP;

int main( int argc, char** argv ) {

	uint16_t targetDimension = 2;
	if( argc >= 2 ) {
		targetDimension = (uint16_t)atoi(argv[1]);
	}

	uint16_t numRBFs = 30;
	if( argc >= 3 ) {
		numRBFs = (uint16_t)atoi(argv[2]);
	}

	// Create a new data set
	DataSystem data;
	//data.maxPoints = 200;

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
	projSecant.targetDimension = targetDimension;
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
	VectorXd parameter(1);
	parameter(0) = 2.5;

	// Compute projected data
	ReducedDataSystem reducedData;
	reducedData.ComputeData(brusselator,data,projSecant.W);

	// Obtain the reduced model
	ModelReducedProducer producer;
	producer.numRBFs = numRBFs;
	ModelReduced reducedModel = producer.ComputeModelReduced(reducedData,data.parameterDimension,data.parameters);

	cout << "Total Cost = " << producer.ComputeTotalCost(reducedModel,reducedData,data.parameters) << endl;

	reducedModel.OutputText("output/reduced.csv");
	//reducedData.WritePointsText("output/p2.5-points.csv");
	//reducedData.WriteVectorsText("output/p2.5-vectors.csv");
	projSecant.WriteBinary("output/projection.bin");
	projSecant.WriteText("output/projection.csv");

	

	return 0;
}
