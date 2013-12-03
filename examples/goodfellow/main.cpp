#include <iostream>
#include <DRDSP/data/data_set.h>
#include <DRDSP/data/secants.h>
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/model_rbf_producer.h>
#include <DRDSP/dynamics/model_reduced_producer.h>
#include <DRDSP/dynamics/generate_data.h>
#include <DRDSP/data/histogram.h>

#include "goodfellow.h"

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

	// The example
	Goodfellow goodfellow(100);
	
	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator dataGenerator(goodfellow);
	dataGenerator.pMin = 4.0;
	dataGenerator.pMax = 5.0;
	dataGenerator.pDelta = (dataGenerator.pMax - dataGenerator.pMin)/10;
	dataGenerator.initial.setRandom();
	dataGenerator.initial -= 0.5 * VectorXd::Ones(goodfellow.dimension);
	dataGenerator.initial *= 2.0;
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 1000;
	dataGenerator.print = 500;
	dataGenerator.rk.dtMax = 0.001;

	DataSystem data = dataGenerator.GenerateDataSystem();
	data.WriteDataSetsCSV("output/orig",".csv");

	// Pre-compute secants
	cout << "Computing secants..." << endl;
	SecantsPreComputed* secants = new SecantsPreComputed [data.numParameters];
	
	for(uint16_t i=0;i<data.numParameters;i++)
		secants[i].ComputeFromData( data.dataSets[i] );

	// Secant culling
	cout << "Culling secants..." << endl;
	SecantsPreComputed* newSecants = new SecantsPreComputed [data.numParameters];
	for(uint16_t i=0;i<data.numParameters;i++)
		newSecants[i] = secants[i].CullSecantsDegrees( 10.0 );

	delete[] secants;

	// Find a projection
	cout << "Finding projection..." << endl;
	ProjSecant projSecant;
	projSecant.targetDimension = options.targetDimension;
	projSecant.targetMinProjectedLength = 0.7;

	// Compute initial condition
	projSecant.GetInitial(data);

	// Optimize over Grassmannian
	projSecant.Find(newSecants,data.numParameters);

	// Print some statistics
	projSecant.AnalyseSecants(newSecants,data.numParameters);

	projSecant.WriteBinary("output/projection.bin");
	projSecant.WriteCSV("output/projection.csv");

	delete[] newSecants;

	// Compute projected data
	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData(goodfellow,data,projSecant.W);
	reducedData.WritePointsCSV("output/p","-points.csv");
	reducedData.WriteVectorsCSV("output/p","-vectors.csv");

	// Obtain the reduced model
	cout << "Computing Reduced Model..." << endl;
	
	ModelReducedProducer producer(options.numRBFs);
	ModelReduced reducedModel = producer.BruteForce(reducedData,data.parameterDimension,data.parameters,options.numIterations);
	
	cout << "Total Cost = " << producer.ComputeTotalCost(reducedModel,reducedData,data.parameters) << endl;
	
	reducedModel.WriteCSV("output/reduced.csv");
	
	//ModelRBFProducer rbfProducer(options.numRBFs);
	//ModelRBF rbfModel = rbfProducer.BruteForce(reducedData.reducedData[0],options.numIterations);

	//cout << "Total Cost = " << rbfProducer.ComputeTotalCost(rbfModel,reducedData.reducedData[0]) << endl;	
	
	// Generate the data
	cout << "Generating Reduced data..." << endl;
	DataGenerator rdataGenerator(reducedModel);
	rdataGenerator.MatchSettings(dataGenerator);
	rdataGenerator.initial = reducedData.reducedData[0].points[0];
	rdataGenerator.tStart = 0.0;

	DataSystem rdata = rdataGenerator.GenerateDataSystem();
	rdata.WriteDataSetsCSV("output/rdata",".csv");

	system("PAUSE");
	return 0;
}
