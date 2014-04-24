#include <iostream>
#include <DRDSP/data/data_set.h>
#include <DRDSP/data/secants.h>
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/model_rbf_producer.h>
#include <DRDSP/dynamics/model_reduced_producer.h>
#include <DRDSP/dynamics/generate_data.h>

#include "kuramoto.h"

using namespace std;
using namespace DRDSP;

struct Options {
	Options() : numIterations(5000), maxPoints(0), targetDimension(4), numRBFs(50) {}
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

	// The kuramoto example
	KuramotoBFlat kuramoto(100);
	
	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator dataGenerator(kuramoto.model);
	dataGenerator.pMin = 2.0;
	dataGenerator.pMax = 2.01;
	dataGenerator.pDelta = 0.002;
	dataGenerator.initial.setRandom();
	dataGenerator.tStart = 50;
	dataGenerator.tInterval = 7;
	dataGenerator.print = 100;
	dataGenerator.rk.dtMax = 0.001;

	DataSystem data = dataGenerator.GenerateDataSystem();
			
	// Embed the data
	cout << "Embedding data..." << endl;
	DataSystem dataEmbedded = kuramoto.embedding.EmbedData(data);

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
	// For this particular example, we use a custom initial condition
	projSecant.W.setZero(kuramoto.embedding.eDim,4);
	
	for(uint32_t i=0;i<(kuramoto.embedding.eDim-2)/2;i++) {
		projSecant.W(2*i,0) = 1;
		projSecant.W(2*i+1,1) = 1;
	}
	projSecant.W(kuramoto.embedding.eDim-2,2) = 1;
	projSecant.W(kuramoto.embedding.eDim-1,3) = 1;
	
	projSecant.W.col(0).normalize();
	projSecant.W.col(1).normalize();

	// Optimize over Grassmannian
	projSecant.Find(newSecants,dataEmbedded.numParameters);

	// Print some statistics
	projSecant.AnalyseSecants(newSecants,dataEmbedded.numParameters);

	projSecant.WriteBinary("output/projection.bin");
	projSecant.WriteCSV("output/projection.csv");

	delete[] newSecants;

	// Compute projected data
	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData(kuramoto,data,projSecant.W);
	reducedData.WritePointsCSV("output/p","-points.csv");
	reducedData.WriteVectorsCSV("output/p","-points.csv");

	// Obtain the reduced model
	cout << "Computing Reduced Model..." << endl;
	
	ModelReducedProducer producer(options.numRBFs);
	ModelReduced reducedModel = producer.BruteForce(reducedData,data.parameterDimension,data.parameters.data(),options.numIterations);
	
	cout << "Total Cost = " << producer.ComputeTotalCost(reducedModel,reducedData,data.parameters.data()) << endl;
	
	reducedModel.WriteCSV("output/reduced.csv");
	
	//ModelRBFProducer rbfProducer(options.numRBFs);
	//ModelRBF rbfModel = rbfProducer.BruteForce(reducedData.reducedData[0],options.numIterations);

	//cout << "Total Cost = " << rbfProducer.ComputeTotalCost(rbfModel,reducedData.reducedData[0]) << endl;	
	
	// Generate the data
	cout << "Generating Reduced data..." << endl;
	DataGenerator rdataGenerator(reducedModel);
	rdataGenerator.MatchSettings(dataGenerator);
	rdataGenerator.initial = reducedData.reducedData[0].points[0];

	DataSystem rdata = rdataGenerator.GenerateDataSystem();
	rdata.WriteDataSetsCSV("output/rdata",".csv");

	system("PAUSE");
	return 0;
}
