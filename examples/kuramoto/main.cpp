#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/rbf_family_producer.h>
#include <DRDSP/dynamics/data_generator.h>
#include <DRDSP/misc.h>
#include "kuramoto.h"

using namespace std;
using namespace DRDSP;

struct Options {
	uint32_t targetDimension, numRBFs, numIterations, numThreads;

	Options() : targetDimension(4), numRBFs(20), numIterations(500), numThreads(4) {}
	
	Options( int argc, char** argv ) : Options() {
		if( argc >= 2 ) targetDimension = (uint32_t)atoi(argv[1]);
		if( argc >= 3 ) numRBFs = (uint32_t)atoi(argv[2]);
		if( argc >= 4 ) numIterations = (uint32_t)atoi(argv[3]);
		if( argc >= 5 ) numThreads = (uint32_t)atoi(argv[4]);
	}
};

typedef Multiquadratic RadialType;

int main( int argc, char** argv ) {

	Options options(argc,argv);
	
	// The kuramoto example
	FamilyEmbedded<KuramotoAFamily,FlatEmbedding> kuramoto(KuramotoAFamily(100),FlatEmbedding(101));
	
	auto parameters = ParameterList( 1.0, 1.1, 6 );

	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator<KuramotoAFamily,KuramotoASolver> dataGenerator(kuramoto.family);
	dataGenerator.initial.setRandom();
	dataGenerator.tStart = 50;
	dataGenerator.tInterval = 7;
	dataGenerator.print = 100;
	dataGenerator.dtMax = 0.001;

	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );
	
	// Embed the data
	cout << "Embedding data..." << endl;
	DataSystem dataEmbedded = EmbedData( kuramoto.embedding, data, options.numThreads );

	// Pre-compute secants
	cout << "Computing secants..." << endl;
	vector<SecantsPreComputed> secants = ComputeSecants( dataEmbedded, options.numThreads );

	// Secant culling
	cout << "Culling secants..." << endl;
	vector<SecantsPreComputed> newSecants = CullSecants( secants, 10.0, options.numThreads );

	secants = vector<SecantsPreComputed>();

	// Find a projection
	cout << "Finding projection..." << endl;
	ProjSecant projSecant( options.targetDimension );

	// Compute initial condition
	// For this particular example, we use a custom initial condition
	projSecant.W.setZero(kuramoto.embedding.eDim,4);
	
	for(uint32_t i=0;i<(kuramoto.embedding.eDim-2)/2;++i) {
		projSecant.W(2*i,0) = 1;
		projSecant.W(2*i+1,1) = 1;
	}
	projSecant.W(kuramoto.embedding.eDim-2,2) = 1;
	projSecant.W(kuramoto.embedding.eDim-1,3) = 1;
	
	projSecant.W.col(0).normalize();
	projSecant.W.col(1).normalize();

	projSecant.Find( newSecants )             // Optimize over Grassmannian
	          .AnalyseSecants( newSecants )   // Print some statistics
	          .WriteBinary("output/projection.bin")
	          .WriteCSV("output/projection.csv");

	newSecants = vector<SecantsPreComputed>();

	// Compute projected data
	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeDataEmbedded( kuramoto, data, projSecant.W, options.numThreads )
	           .WritePointsCSV("output/p","-points.csv")
	           .WriteVectorsCSV("output/p","-points.csv");

	// Obtain the reduced model
	cout << "Computing Reduced Model..." << endl;
	
	RBFFamilyProducer<RadialType> producer(options.numRBFs);
	auto reducedModel = producer.BruteForce(reducedData,data.parameterDimension,data.parameters,options.numIterations);
	
	cout << "Total Cost = " << producer.ComputeTotalCost(reducedModel,reducedData,data.parameters) << endl;
	
	reducedModel.WriteCSV("output/reduced.csv");
	
	//ModelRBFProducer rbfProducer(options.numRBFs);
	//ModelRBF rbfModel = rbfProducer.BruteForce(reducedData.reducedData[0],options.numIterations);

	//cout << "Total Cost = " << rbfProducer.ComputeTotalCost(rbfModel,reducedData.reducedData[0]) << endl;	
	
	// Generate the data
	cout << "Generating Reduced data..." << endl;
	DataGenerator<RBFFamily<RadialType>> rdataGenerator(reducedModel);
	rdataGenerator.MatchSettings(dataGenerator);
	rdataGenerator.tStart = 0.0;

	DataSystem rdata = rdataGenerator.GenerateUsingInitials( parameters, reducedData, options.numThreads );
	rdata.WriteDataSetsCSV("output/rdata",".csv");

	Compare( reducedData, rdata );

	cout << "Press any key to continue . . . "; cin.get();
}
