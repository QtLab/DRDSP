#include <iostream>
#include <DRDSP/data/data_set.h>
#include <DRDSP/data/secants.h>
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/model_reduced_producer.h>
#include <DRDSP/dynamics/generate_data.h>
#include <DRDSP/dynamics/bifurcation.h>
#include <DRDSP/dynamics/rk.h>

#include "rossler.h"

using namespace std;
using namespace DRDSP;

struct Options {
	Options() : numIterations(1000), maxPoints(0), targetDimension(3), numRBFs(40) {}
	uint32_t numIterations, maxPoints;
	uint32_t targetDimension, numRBFs;
};

Options GetOptions( int argc, char** argv ) {
	Options options;

	if( argc >= 2 ) options.targetDimension = (uint16_t)atoi(argv[1]);
	if( argc >= 3 ) options.numRBFs = (uint16_t)atoi(argv[2]);
	if( argc >= 4 ) options.maxPoints = (uint32_t)atoi(argv[3]);
	if( argc >= 5 ) options.numIterations = (uint32_t)atoi(argv[4]);

	return options;
}

void Compare( const ReducedDataSystem& reducedData, const DataSystem& rdata ) {

	ofstream out("output/comparison.csv");
	out << "Parameter,RMS,Max,Differences" << endl;
	for(uint32_t i=0;i<reducedData.numParameters;++i) {
		DataComparisonResult r = CompareData( reducedData.reducedData[i].points, rdata.dataSets[i].points );
		cout << "Parameter " << rdata.parameters[i] << endl;
		cout << "RMS: " << r.rmsDifference << endl;
		cout << "Max: " << r.maxDifference << endl;
		
		out << rdata.parameters[i] << ",";
		out << r.rmsDifference << ",";
		out << r.maxDifference << ",";
		for( const auto& x : r.differences )
			out << x << ",";
		out << endl;
	}

}

int main( int argc, char** argv ) {

	Options options = GetOptions(argc,argv);

	// The example
	RosslerFamily rossler;

	auto parameters = ParameterList( 4.0, 8.8, 21 );

	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator<RosslerFamily> dataGenerator;
	dataGenerator.initial.setRandom();
	dataGenerator.initial -= 0.5 * VectorXd::Ones(rossler.dimension);
	dataGenerator.initial *= 2.0;
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 100;
	dataGenerator.dtMax = 0.0005;
	dataGenerator.print = 1000;
	
	DataSystem data = dataGenerator.GenerateDataSystem( parameters );
	data.WriteDataSetsCSV("output/orig",".csv");

	// Compute projected data
	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData( rossler, data, MatrixXd::Identity(3,3) );
	reducedData.WritePointsCSV( "output/p", "-points.csv" );
	reducedData.WriteVectorsCSV( "output/p", "-vectors.csv" );

	// Obtain the reduced model
	cout << "Computing Reduced Model..." << endl;
	
	ModelReducedProducer<> producer(options.numRBFs);
	auto reducedModel = producer.BruteForce(reducedData,data.parameterDimension,data.parameters.data(),options.numIterations);
	
	cout << "Total Cost = " << producer.ComputeTotalCost(reducedModel,reducedData,data.parameters.data()) << endl;
	
	reducedModel.WriteCSV("output/reduced.csv");

	// Generate the data
	cout << "Generating Reduced data..." << endl;
	DataGenerator<ModelReduced<>> rdataGenerator(reducedModel);
	rdataGenerator.MatchSettings(dataGenerator);
	rdataGenerator.initial = reducedData.reducedData[0].points[0];
	rdataGenerator.tStart = 0.0;

	DataSystem rdata = rdataGenerator.GenerateUsingInitials( parameters, reducedData );
	rdata.WriteDataSetsCSV("output/rdata",".csv");

	Compare( reducedData, rdata );

	// Bifurcation diagrams
	cout << "Generating Bifurcation Diagram..." << endl;
	{
		BifurcationDiagramGenerator<RosslerFamily> bifurcationGenerator;
		bifurcationGenerator.startTime = 1000.0;
		bifurcationGenerator.endTime = 2000.0;
		bifurcationGenerator.dt = 0.001;
		bifurcationGenerator.minParam = parameters[0];
		bifurcationGenerator.maxParam = parameters[20];
		bifurcationGenerator.numParams = 1280;
		bifurcationGenerator.initial = dataGenerator.initial;
		bifurcationGenerator.Generate( rossler,
			[]( const VectorXd& x, const VectorXd& y ) { return x[1] > 0.0 && y[1] < 0.0; },
			[]( const VectorXd& x, const VectorXd& y ) { return (x[0]+y[0])*0.5; }
		).WriteBitmap("output/bifurcation-orig.bmp",720);
	}
	{
		BifurcationDiagramGenerator<ModelReduced<>> bifurcationGenerator;
		bifurcationGenerator.startTime = 1000.0;
		bifurcationGenerator.endTime = 2000.0;
		bifurcationGenerator.dt = 0.001;
		bifurcationGenerator.minParam = parameters[0];
		bifurcationGenerator.maxParam = parameters[20];
		bifurcationGenerator.numParams = 1280;
		bifurcationGenerator.initial = dataGenerator.initial;
		bifurcationGenerator.Generate( reducedModel,
			[]( const VectorXd& x, const VectorXd& y ) { return x[1] > 0.0 && y[1] < 0.0; },
			[]( const VectorXd& x, const VectorXd& y ) { return (x[0]+y[0])*0.5; }
		).WriteBitmap("output/bifurcation-red.bmp",720);
	}

	system("PAUSE");
	return 0;
}
