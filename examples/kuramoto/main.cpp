#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/rbf_family_producer.h>
#include <DRDSP/dynamics/data_generator.h>
#include <DRDSP/misc.h>
#include "kuramoto.h"

using namespace std;
using namespace DRDSP;

struct Options {
	uint32_t targetDimension = 4,
	         numRBFs = 20,
	         numIterations = 500,
	         numThreads = 4;
	
	Options( int argc, char** argv ) {
		if( argc >= 2 ) targetDimension = (uint32_t)atoi(argv[1]);
		if( argc >= 3 )         numRBFs = (uint32_t)atoi(argv[2]);
		if( argc >= 4 )   numIterations = (uint32_t)atoi(argv[3]);
		if( argc >= 5 )      numThreads = (uint32_t)atoi(argv[4]);
	}
};

typedef RBF<Multiquadratic> RBFType;

int main( int argc, char** argv ) {

	Options options(argc,argv);
	
	FamilyEmbedded<KuramotoAFamily,FlatEmbedding> kuramoto(KuramotoAFamily(100),FlatEmbedding(101));
	
	auto parameters = ParameterList( 1.0, 1.1, 6 );

	cout << "Generating data..." << endl;
	DataGenerator<KuramotoAFamily,KuramotoASolver> dataGenerator( kuramoto.family );
	dataGenerator.initial.setRandom();
	dataGenerator.tStart = 50;
	dataGenerator.tInterval = 7;
	dataGenerator.print = 100;
	dataGenerator.dtMax = 0.001;

	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );
	
	cout << "Embedding data..." << endl;
	DataSystem dataEmbedded = EmbedData( kuramoto.embedding, data, options.numThreads );

	cout << "Computing secants..." << endl;
	vector<Secants> secants = ComputeSecants( dataEmbedded, 10.0, options.numThreads );

	cout << "Finding projection..." << endl;
	ProjSecant projSecant( options.targetDimension );

	// For this particular example, we use a custom initial condition
	projSecant.W.setZero(kuramoto.embedding.embedDim,4);
	
	for(uint32_t i=0;i<(kuramoto.embedding.embedDim-2)/2;++i) {
		projSecant.W(2*i,0) = 1;
		projSecant.W(2*i+1,1) = 1;
	}
	projSecant.W(kuramoto.embedding.embedDim-2,2) = 1;
	projSecant.W(kuramoto.embedding.embedDim-1,3) = 1;
	
	projSecant.W.col(0).normalize();
	projSecant.W.col(1).normalize();

	projSecant.Find( secants )             // Optimize over Grassmannian
	          .AnalyseSecants( secants )   // Print some statistics
	          .WriteBinary("output/projection.bin")
	          .WriteCSV("output/projection.csv");

	secants = vector<Secants>();

	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeDataEmbedded( kuramoto, data, projSecant.W, options.numThreads )
	           .WritePointsCSV("output/p","-points.csv")
	           .WriteVectorsCSV("output/p","-points.csv");

	cout << "Computing Reduced Family..." << endl;
	
	RBFFamilyProducer<RBFType> producer( options.numRBFs );
	auto reducedFamily = producer.BruteForce( reducedData,
	                                          data.parameters,
	                                          options.numIterations,
	                                          options.numThreads );
	
	cout << "Total Cost = "
	     << producer.ComputeTotalCost( reducedFamily, reducedData, data.parameters )
	     << endl;
	
	//reducedFamily.WriteCSV("output/reduced.csv");
	
	//ModelRBFProducer rbfProducer( options.numRBFs );
	//ModelRBF rbfModel = rbfProducer.BruteForce( reducedData.reducedData[0], options.numIterations );

	//cout << "Total Cost = " << rbfProducer.ComputeTotalCost( rbfModel, reducedData.reducedData[0] ) << endl;	

	cout << "Generating Reduced data..." << endl;
	auto rdataGenerator = MakeDataGenerator( reducedFamily );
	rdataGenerator.MatchSettings( dataGenerator );
	rdataGenerator.tStart = 0.0;

	DataSystem rdata = rdataGenerator.GenerateUsingInitials( parameters, reducedData, options.numThreads );
	rdata.WriteCSV("output/rdata",".csv");

	Compare( reducedData, rdata );

	cout << "Press enter to continue . . . "; cin.get();
}
