#pragma warning ( disable : 4503 )
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/rbf_family_producer.h>
#include <DRDSP/dynamics/data_generator.h>
#include "dynamo.h"

using namespace std;
using namespace DRDSP;

struct Options {
	uint32_t targetDimension, numRBFs, numIterations, numThreads;

	Options() : targetDimension(2), numRBFs(30), numIterations(500), numThreads(3) {}
	
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

	// The example
	DynamoFamily dynamo;

	auto parameters = ParameterList( 1.2, 2.0, 3 );
	
	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator<DynamoFamily> dataGenerator(dynamo);
	dataGenerator.initial.setRandom(dynamo.dimension);
	dataGenerator.tStart = 0;
	dataGenerator.tInterval = 0.00012;
	dataGenerator.print = 200;
	dataGenerator.dtMax = 1.0e-7;

	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );
	data.WriteDataSetsBinary("output/orig",".bin");

	// Pre-compute secants
	cout << "Computing secants..." << endl;
	vector<SecantsPreComputed> secants = ComputeSecants( data, options.numThreads );

	// Secant culling
	cout << "Culling secants..." << endl;
	vector<SecantsPreComputed> newSecants = CullSecants( secants, 10.0, options.numThreads );

	secants = vector<SecantsPreComputed>();

	// Find a projection
	cout << "Finding projection..." << endl;
	ProjSecant projSecant;
	projSecant.targetDimension = options.targetDimension;
	projSecant.targetMinProjectedLength = 0.7;

	// Compute initial condition
	projSecant.GetInitial(data);

	// Optimize over Grassmannian
	projSecant.Find( newSecants );

	// Print some statistics
	projSecant.AnalyseSecants( newSecants );

	projSecant.WriteBinary("output/projection.bin");
	projSecant.WriteCSV("output/projection.csv");

	newSecants = vector<SecantsPreComputed>();

	// Compute projected data
	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData( dynamo, data, projSecant.W, options.numThreads );
	reducedData.WritePointsCSV("output/p","-points.csv");
	reducedData.WriteVectorsCSV("output/p","-vectors.csv");

	// Obtain the reduced model
	cout << "Computing Reduced Model..." << endl;
	
	RBFFamilyProducer<RadialType> producer(options.numRBFs);
	auto reducedModel = producer.BruteForce(reducedData,data.parameterDimension,data.parameters,options.numIterations);
	
	cout << "Total Cost = " << producer.ComputeTotalCost(reducedModel,reducedData,data.parameters) << endl;
	
	reducedModel.WriteCSV("output/reduced.csv");

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

