#include <iostream>
#include <DRDSP/data/data_set.h>
#include <DRDSP/data/secants.h>
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/model_reduced_producer.h>
#include <DRDSP/dynamics/generate_data.h>
#include "brusselator.h"

using namespace std;
using namespace DRDSP;

struct Options {
	Options() : numIterations(1000), maxPoints(0), targetDimension(2), numRBFs(30) {}
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

	// Dynamics
	BrusselatorFamily brusselator;

	auto parameters = ParameterList( 2.5, 2.6, 5 );

	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator<BrusselatorFamily> dataGenerator;
	dataGenerator.initial.setRandom(brusselator.dimension);
	dataGenerator.tStart = 30;
	dataGenerator.tInterval = 6.8;
	dataGenerator.print = 100;
	dataGenerator.dtMax = 0.001;

	DataSystem data = dataGenerator.GenerateDataSystem( parameters );

	// Pre-compute secants
	SecantsPreComputed* secants = new SecantsPreComputed [data.numParameters];
	SecantsPreComputed* newSecants = new SecantsPreComputed [data.numParameters];
	for(uint16_t i=0;i<data.numParameters;i++)
		secants[i].ComputeFromData( data.dataSets[i] );

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

	projSecant.WriteBinary("output/projection.bin");
	projSecant.WriteCSV("output/projection.csv");

	
	// Compute projected data
	cout << endl << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData( brusselator, data, projSecant.W );

	reducedData.WritePointsCSV("output/p","-points.csv");
	reducedData.WriteVectorsCSV("output/p","-vectors.csv");


	// Obtain the reduced model
	cout << endl << "Computing Reduced Model..." << endl;
	ModelReducedProducer<RadialType> producer(options.numRBFs);
	auto reducedModel = producer.BruteForce(reducedData,data.parameterDimension,data.parameters.data(),options.numIterations);

	cout << "Total Cost = " << producer.ComputeTotalCost(reducedModel,reducedData,data.parameters.data()) << endl;

	reducedModel.WriteCSV("output/reduced.csv");
	

	// Generate the data
	cout << "Simulating the reduced model..." << endl;
	DataGenerator<ModelReduced<RadialType>> rdataGenerator(reducedModel);
	rdataGenerator.MatchSettings(dataGenerator);
	rdataGenerator.tStart = 0.0;
	DataSystem rdata = rdataGenerator.GenerateUsingInitials( parameters, reducedData );
	rdata.WriteDataSetsCSV("output/rdata",".csv");

	Compare( reducedData, rdata );

	system("PAUSE");
}
