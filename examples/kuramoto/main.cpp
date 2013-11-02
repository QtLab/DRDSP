#include <iostream>
#include <DRDSP/data/data_set.h>
#include <DRDSP/data/secants.h>
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/model_reduced_producer.h>
#include <DRDSP/dynamics/generate_data.h>

#include "kuramoto.h"

using namespace std;
using namespace DRDSP;

struct Options {
	Options() : numIterations(1000), maxPoints(0), targetDimension(10), numRBFs(50) {}
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

	// The ks example
	KuramotoFlat kuramoto(100);

	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator dataGenerator(kuramoto.model);
	dataGenerator.pMin = 1.0;
	dataGenerator.pMax = 2.0;
	dataGenerator.pDelta = 0.1;
	dataGenerator.initial.setRandom();
	dataGenerator.tStart = 50;
	dataGenerator.tInterval = 7;
	dataGenerator.print = 100;
	dataGenerator.rk.dtMax = 0.001;

	DataSystem data = dataGenerator.GenerateDataSystem();
	
	data.dataSets[0].WriteCSV("output/orig.csv");
		
	// Embed the data
	cout << "Embedding data..." << endl;
	DataSystem dataEmbedded = kuramoto.embedding.EmbedData(data);

	dataEmbedded.dataSets[0].WriteCSV("output/embed.csv");

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
	reducedData.ComputeData(kuramoto,data,projSecant.W);

	// Obtain the reduced model
	cout << "Computing Reduced Model..." << endl;
	ModelReducedProducer producer;
	producer.numRBFs = options.numRBFs;
	ModelReduced reducedModel = producer.BruteForce(reducedData,data.parameterDimension,data.parameters,options.numIterations);

	cout << "Total Cost = " << producer.ComputeTotalCost(reducedModel,reducedData,data.parameters) << endl;

	reducedModel.WriteCSV("output/reduced.csv");
	reducedData.reducedData[0].WritePointsText("output/p1-points.csv");
	projSecant.WriteBinary("output/projection.bin");
	projSecant.WriteCSV("output/projection.csv");

	// Generate the data
	cout << "Generating Reduced data..." << endl;
	DataGenerator rdataGenerator(reducedModel);
	rdataGenerator.MatchSettings(dataGenerator);
	rdataGenerator.initial = reducedData.reducedData[0].points[0];

	DataSystem rdata = rdataGenerator.GenerateDataSystem();

	rdata.dataSets[0].WriteCSV("output/rdata.csv");

	system("PAUSE");
	return 0;
}
