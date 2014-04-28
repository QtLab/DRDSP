#include <iostream>
#include <DRDSP/data/data_set.h>
#include <DRDSP/data/secants.h>
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/model_reduced_producer.h>
#include <DRDSP/dynamics/generate_data.h>

#include "pendulum.h"

using namespace std;
using namespace DRDSP;

struct Options {
	Options() : numIterations(1000), maxPoints(0), targetDimension(4), numRBFs(100) {}
	uint32_t numIterations, maxPoints, targetDimension, numRBFs;
};

Options GetOptions( int argc, char** argv ) {
	Options options;

	if( argc >= 2 ) options.targetDimension = (uint32_t)atoi(argv[1]);
	if( argc >= 3 ) options.numRBFs = (uint32_t)atoi(argv[2]);
	if( argc >= 4 ) options.maxPoints = (uint32_t)atoi(argv[3]);
	if( argc >= 5 ) options.numIterations = (uint32_t)atoi(argv[4]);

	return options;
}

void Compare( const ReducedDataSystem& reducedData, const DataSystem& rdata ) {

	ofstream out("output/comparison.csv");
	out << "Parameter,RMS,Max,MaxMin,Differences" << endl;
	for(uint32_t i=0;i<reducedData.numParameters;++i) {
		DataComparisonResult r = CompareData( reducedData.reducedData[i].points, rdata.dataSets[i].points );
		cout << "Parameter " << rdata.parameters[i] << endl;
		cout << "RMS: " << r.rmsDifference << endl;
		cout << "Max: " << r.maxDifference << endl;
		cout << "MaxMin: " << r.maxMinDifference << endl;

		out << rdata.parameters[i] << ",";
		out << r.rmsDifference << ",";
		out << r.maxDifference << ",";
		out << r.maxMinDifference << ",";
		for( const auto& x : r.differences )
			out << x << ",";
		out << endl;
	}

}

typedef Multiquadratic RadialType;

int main( int argc, char** argv ) {

	Options options = GetOptions(argc,argv);

	// The pendulum example
	FamilyEmbedded<PendulumFamily,FlatEmbedding> pendulum;

	auto parameters = ParameterList( 1.8, 1.82, 5 );

	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator<PendulumFamily,PendulumSolver> dataGenerator;
	dataGenerator.initial.setZero(pendulum.family.dimension);
	dataGenerator.initial(1) = -0.422;
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 7;
	dataGenerator.print = 100;
	dataGenerator.dtMax = 0.001;

	DataSystem data = dataGenerator.GenerateDataSystem( parameters );
	data.WriteDataSetsCSV("output/orig",".csv");
	
	// Embed the data
	cout << "Embedding data..." << endl;
	DataSystem dataEmbedded = EmbedData(pendulum.embedding,data);
	
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
	projSecant.Find( newSecants, dataEmbedded.numParameters );

	// Print some statistics
	projSecant.AnalyseSecants( newSecants, dataEmbedded.numParameters );

	delete[] newSecants;

	projSecant.WriteBinary("output/projection.bin");
	projSecant.WriteCSV("output/projection.csv");


	// Compute projected data
	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeDataEmbedded( pendulum, data, projSecant.W );

	reducedData.WritePointsCSV("output/p","-points.csv");
	reducedData.WriteVectorsCSV("output/p","-vectors.csv");

	
	// Obtain the reduced model
	cout << "Computing Reduced Model..." << endl;
	ModelReducedProducer<RadialType> producer(options.numRBFs);
	auto reducedModel = producer.BruteForce( reducedData,
	                                         data.parameterDimension,
	                                         data.parameters.data(),
	                                         options.numIterations );

	cout << "Total Cost = "
		 << producer.ComputeTotalCost( reducedModel, reducedData, data.parameters.data() )
		 << endl;

	reducedModel.WriteCSV("output/reduced.csv");

	// Generate the data
	cout << "Simulating the reduced model..." << endl;
	DataGenerator<ModelReduced<RadialType>> rdataGenerator(reducedModel);
	rdataGenerator.MatchSettings(dataGenerator);
	rdataGenerator.tStart = 0.0;
	DataSystem rdata = rdataGenerator.GenerateUsingInitials( parameters, reducedData );

	Compare( reducedData, rdata );

	system("PAUSE");
}
