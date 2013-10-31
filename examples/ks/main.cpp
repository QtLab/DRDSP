#include <iostream>
#include <DRDSP/data/data_set.h>
#include <DRDSP/data/secants.h>
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/model_reduced_producer.h>
#include <DRDSP/dynamics/generate_data.h>

#include "ks.h"

using namespace std;
using namespace DRDSP;

struct Options {
	Options() : numIterations(1000), maxPoints(0), targetDimension(4), numRBFs(50) {}
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

	// The pendulum example
	KSFlat ks;

	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator dataGenerator(ks.model);
	dataGenerator.pMin = 1.8;
	dataGenerator.pMax = 1.8201;
	dataGenerator.pDelta = 0.005;
	dataGenerator.initial(0) = 0.3;
	dataGenerator.initial(1) = 0.3;
	dataGenerator.tStart = 10000;
	dataGenerator.tInterval = 3.5;
	dataGenerator.print = 100;
	dataGenerator.rk.dtMax = 0.001;

	DataSystem data = dataGenerator.GenerateDataSystem();
	
	// Embed the data
	cout << "Embedding data..." << endl;
	DataSystem dataEmbedded = ks.embedding.EmbedData(data);
	
	// Pre-compute secants
	cout << "Computing secants..." << endl;
	SecantsPreComputed* secants = new SecantsPreComputed [dataEmbedded.numParameters];
	
	for(uint16_t i=0;i<dataEmbedded.numParameters;i++)
		secants[i].ComputeFromData( dataEmbedded.dataSets[i] );

	// Secant culling
	cout << "Culling secants..." << endl;
	SecantsPreComputed* newSecants = new SecantsPreComputed [dataEmbedded.numParameters];
	for(uint16_t i=0;i<dataEmbedded.numParameters;i++)
		newSecants[i] = secants[i].CullSecantsDegrees( 10.0 );

	delete[] secants;

	// Find a projection
	cout << "Finding projection..." << endl;
	ProjSecant projSecant;
	projSecant.targetDimension = options.targetDimension;
	projSecant.targetMinProjectedLength = 0.7;

	// Compute initial condition
	projSecant.GetInitial(dataEmbedded);

	// Optimize over Grassmannian
	projSecant.Find(newSecants,dataEmbedded.numParameters);

	// Print some statistics
	projSecant.AnalyseSecants(newSecants,dataEmbedded.numParameters);

	delete[] newSecants;

	// Compute projected data
	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData(ks,data,projSecant.W);

	// Obtain the reduced model
	cout << "Computing Reduced Model..." << endl;
	ModelReducedProducer producer;
	producer.numRBFs = options.numRBFs;
	ModelReduced reducedModel = producer.BruteForce(reducedData,data.parameterDimension,data.parameters,options.numIterations);

	cout << "Total Cost = " << producer.ComputeTotalCost(reducedModel,reducedData,data.parameters) << endl;

	reducedModel.OutputText("output/reduced.csv");
	//reducedData.WritePointsText("output/p2.5-points.csv");
	//reducedData.WriteVectorsText("output/p2.5-vectors.csv");
	projSecant.WriteBinary("output/projection.bin");
	projSecant.WriteText("output/projection.csv");

	// Generate the data
	cout << "Generating Reduced data..." << endl;
	DataGenerator rdataGenerator(reducedModel);
	rdataGenerator.pMin = 1.8;
	rdataGenerator.pMax = 1.8251;
	rdataGenerator.pDelta = 0.005;
	rdataGenerator.initial = reducedData.reducedData[0].points[0];
	rdataGenerator.tStart = 0;
	rdataGenerator.tInterval = 3.5;
	rdataGenerator.print = 100;
	rdataGenerator.rk.dtMax = 0.001;

	DataSystem rdata = rdataGenerator.GenerateDataSystem();

	rdata.dataSets[0].WriteText("output/p1.8.csv");

	system("PAUSE");
	return 0;
}
