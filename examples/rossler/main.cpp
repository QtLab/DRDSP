#pragma warning( disable : 4530 )
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/rbf_family_producer.h>
#include <DRDSP/dynamics/data_generator.h>
#include <DRDSP/dynamics/bifurcation.h>
#include <DRDSP/dynamics/monodromy.h>
#include "rossler.h"

using namespace std;
using namespace DRDSP;

typedef PolyharmonicSpline3 RadialType;

struct Options {
	uint32_t targetDimension, numRBFs, numIterations, numThreads;

	Options() : targetDimension(3), numRBFs(40), numIterations(150), numThreads(4) {}
	
	Options( int argc, char** argv ) : Options() {
		if( argc >= 2 ) targetDimension = (uint32_t)atoi(argv[1]);
		if( argc >= 3 ) numRBFs = (uint32_t)atoi(argv[2]);
		if( argc >= 4 ) numIterations = (uint32_t)atoi(argv[3]);
		if( argc >= 5 ) numThreads = (uint32_t)atoi(argv[4]);
	}
};

void ComputeReduced( const Options& );
void SimulateReduced( const Options& );
void OriginalFloquet( const Options& );
void ReducedFloquet( const Options& );

int main( int argc, char** argv ) {
	Options options(argc,argv);

	//ComputeReduced( options );
	SimulateReduced( options );
	//OriginalFloquet( options );
	//ReducedFloquet( options );

	cout << "Press any key to continue . . . "; cin.get();
}

void ReducedFloquet( const Options& options ) {

	auto parameters = ParameterList( 4.0, 8.8, 21 );

	vector<double> periods = {
		6.0, 6.0, 6.0, 6.01, 6.02, 6.03,
		12.05, 12.05, 12.06, 12.07, 12.08, 12.08, 12.09, 12.10, 12.10, 12.11,
		24.23, 24.24, 24.25, 48.51,
		100.0
	};

	RBFFamily<RadialType> reducedModel;
	reducedModel.ReadText("reduced.txt");

	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator<RBFFamily<RadialType>> dataGenerator(reducedModel);
	dataGenerator.initial = Vector3d(7.0,0.0,0.5);
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 100;
	dataGenerator.dtMax = 0.001;
	dataGenerator.print = 10000;
	
	DataSystem data = dataGenerator.GenerateDataSystem( parameters, periods, options.numThreads );
	data.WriteDataSetsCSV("output/orig",".csv");
	
	cout << "Computing Floquet multipliers..." << endl;
	ofstream out("output/floquet.csv");
	for(int i=0;i<parameters.size();++i) {
		double dt = periods[i] / (dataGenerator.print-1);
		auto floquet = ComputeFloquetMultipliers( reducedModel(parameters[i]), data.dataSets[i].points, dt );
		out << parameters[i] << ",";
		for(int j=0;j<floquet.size();++j)
			out << floquet[j].real() << "," << floquet[j].imag() << ",";
		out << endl;
	}
}


void OriginalFloquet( const Options& options ) {

	// The example
	RosslerFamily rossler;

	auto parameters = ParameterList( 4.0, 8.8, 21 );

	vector<double> periods = {
		6.0, 6.0, 6.0, 6.01, 6.02, 6.03,
		12.05, 12.05, 12.06, 12.07, 12.08, 12.08, 12.09, 12.10, 12.10, 12.11,
		24.23, 24.24, 24.25, 48.51,
		100.0
	};

	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator<RosslerFamily> dataGenerator;
	dataGenerator.initial = Vector3d(7.0,0.0,0.5);
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 100;
	dataGenerator.dtMax = 0.001;
	dataGenerator.print = 10000;
	
	DataSystem data = dataGenerator.GenerateDataSystem( parameters, periods, options.numThreads );
	data.WriteDataSetsCSV("output/orig",".csv");
	
	cout << "Computing Floquet multipliers..." << endl;
	ofstream out("output/floquet.csv");
	for(int i=0;i<parameters.size();++i) {
		double dt = periods[i] / (dataGenerator.print-1);
		auto floquet = ComputeFloquetMultipliers( rossler(parameters[i]), data.dataSets[i].points, dt );
		out << parameters[i] << ",";
		for(int j=0;j<floquet.size();++j)
			out << floquet[j].real() << "," << floquet[j].imag() << ",";
		out << endl;
	}
}

void ComputeReduced( const Options& options ) {

	// The example
	RosslerFamily rossler;

	auto parameters = ParameterList( 4.0, 8.8, 21 );

	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator<RosslerFamily> dataGenerator;
	dataGenerator.initial = Vector3d(7.0,0.0,0.5);
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 100;
	dataGenerator.dtMax = 0.001;
	dataGenerator.print = 1000;
	
	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );
	data.WriteDataSetsCSV("output/orig",".csv");

	// Compute projected data
	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData( rossler, data, MatrixXd::Identity(3,3), options.numThreads )
	           .WritePointsCSV( "output/p", "-points.csv" )
	           .WriteVectorsCSV( "output/p", "-vectors.csv" )
	           .WriteDerivativesCSV( "output/p", "-derivatives.csv" );

	// Obtain the reduced model
	cout << "Computing Reduced Model..." << endl;
	
	RBFFamilyProducer<RadialType> producer(options.numRBFs);
	producer.boxScale = 1.6;
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

	// Bifurcation diagrams
	cout << "Generating Bifurcation Diagram..." << endl;
	{
		BifurcationDiagramGenerator<RosslerFamily> bifurcationGenerator;
		bifurcationGenerator.tStart = 500.0;
		bifurcationGenerator.tInterval = 500.0;
		bifurcationGenerator.dt = 0.001;
		bifurcationGenerator.dtMax = 0.001;
		bifurcationGenerator.pMin = parameters.front();
		bifurcationGenerator.pMax = parameters.back();
		bifurcationGenerator.pCount = 1280;
		bifurcationGenerator.initial = dataGenerator.initial;
		bifurcationGenerator.Generate( rossler,
			[]( const VectorXd& x, const VectorXd& y ) { return x[1] > 0.0 && y[1] < 0.0; },
			[]( const VectorXd& x, const VectorXd& y ) { return (x[0]+y[0])*0.5; },
			options.numThreads ).WriteBitmap( "output/bifurcation-orig.bmp", 720 );
	}
	{
		BifurcationDiagramGenerator<RBFFamily<RadialType>> bifurcationGenerator;
		bifurcationGenerator.tStart = 500.0;
		bifurcationGenerator.tInterval = 500.0;
		bifurcationGenerator.dt = 0.001;
		bifurcationGenerator.dtMax = 0.001;
		bifurcationGenerator.pMin = parameters.front();
		bifurcationGenerator.pMax = parameters.back();
		bifurcationGenerator.pCount = 1280;
		bifurcationGenerator.initial = dataGenerator.initial;
		bifurcationGenerator.Generate( reducedModel,
			[]( const VectorXd& x, const VectorXd& y ) { return x[1] > 0.0 && y[1] < 0.0; },
			[]( const VectorXd& x, const VectorXd& y ) { return (x[0]+y[0])*0.5; },
			options.numThreads ).WriteBitmap( "output/bifurcation-red.bmp", 720 );
	}
}

void SimulateReduced( const Options& options ) {

	// The example
	RosslerFamily rossler;

	auto parameters = ParameterList( 4.0, 8.8, 21 );

	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator<RosslerFamily> dataGenerator;
	dataGenerator.initial = Vector3d(7.0,0.0,0.5);
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 100;
	dataGenerator.dtMax = 0.001;
	dataGenerator.print = 1000;

	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );
	data.WriteDataSetsCSV("output/orig",".csv");

	// Compute projected data
	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData( rossler, data, MatrixXd::Identity(3,3), options.numThreads );

	RBFFamily<RadialType> reducedModel;
	reducedModel.ReadText("reduced.txt");

	RBFFamilyProducer<RadialType> producer( reducedModel.model.numRBFs );

	cout << "Total Cost = " << producer.ComputeTotalCost(reducedModel,reducedData,data.parameters) << endl;

	producer.Fit( reducedModel, reducedData, rossler.parameterDimension, data.parameters );
	
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

	BifurcationDiagramGenerator<RBFFamily<RadialType>> bifurcationGenerator;
	bifurcationGenerator.tStart = 500.0;
	bifurcationGenerator.tInterval = 500.0;
	bifurcationGenerator.dt = 0.001;
	bifurcationGenerator.dtMax = 0.001;
	bifurcationGenerator.pMin = parameters.front();
	bifurcationGenerator.pMax = parameters.back();
	bifurcationGenerator.pCount = 1280;
	bifurcationGenerator.initial = dataGenerator.initial;
	bifurcationGenerator.Generate( reducedModel,
		[]( const VectorXd& x, const VectorXd& y ) { return x[1] > 0.0 && y[1] < 0.0; },
		[]( const VectorXd& x, const VectorXd& y ) { return (x[0]+y[0])*0.5; },
		options.numThreads ).WriteBitmap( "output/bifurcation-red.bmp", 720 );

}
