#pragma warning( disable : 4530 )
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/rbf_family_producer.h>
#include <DRDSP/dynamics/data_generator.h>
#include <DRDSP/dynamics/bifurcation.h>
#include <DRDSP/dynamics/monodromy.h>
#include <DRDSP/misc.h>
#include "rossler.h"

using namespace std;
using namespace DRDSP;

typedef RBF<PolyharmonicSpline<3>> RBFType;

struct Options {
	uint32_t targetDimension, numRBFs, numIterations, numThreads;

	Options() : targetDimension(3), numRBFs(40), numIterations(250), numThreads(4) {}
	
	Options( int argc, char** argv ) : Options() {
		if( argc >= 2 ) targetDimension = (uint32_t)atoi(argv[1]);
		if( argc >= 3 )         numRBFs = (uint32_t)atoi(argv[2]);
		if( argc >= 4 )   numIterations = (uint32_t)atoi(argv[3]);
		if( argc >= 5 )      numThreads = (uint32_t)atoi(argv[4]);
	}
};

void ComputeReduced( const Options& );
void SimulateReduced( const Options& );
void OriginalFloquet( const Options& );
void ReducedFloquet( const Options& );
void Test( const Options& );
void Test2( const Options& );

int main( int argc, char** argv ) {
	Options options(argc,argv);

	ComputeReduced( options );
	//SimulateReduced( options );
	//OriginalFloquet( options );
	//ReducedFloquet( options );
	//Test( options );
	//Test2( options );

	cout << "Press any key to continue . . . "; cin.get();
}

void ReducedFloquet( const Options& options ) {

	auto parameters = ParameterList( 4.0, 8.8, 21 );

	RBFFamily<RBFType> reducedFamily;
	reducedFamily.ReadText("reduced.txt");

	cout << "Generating data..." << endl;
	DataGenerator<RBFFamily<RBFType>> dataGenerator( reducedFamily );
	dataGenerator.initial = Vector3d(7.0,0.0,0.5);
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 100;
	dataGenerator.dtMax = 0.001;
	dataGenerator.print = 10000;
	
	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );

	vector<double> periods(parameters.size());
	for(size_t i=0;i<parameters.size();++i) {
		periods[i] = DetectPeriod(
			data.dataSets[i].points,
			reducedFamily(parameters[i])(data.dataSets[i][0]),
			dataGenerator.tInterval / (dataGenerator.print-1),
			0.001
		);
		cout << periods[i] << endl;
	}

	data = dataGenerator.GenerateDataSystem( parameters, periods, options.numThreads );
	data.WriteDataSetsCSV("output/orig",".csv");
	
	cout << "Computing Floquet multipliers..." << endl;
	ofstream out("output/floquet.csv");
	for(size_t i=0;i<parameters.size();++i) {
		double dt = periods[i] / (dataGenerator.print-1);
		auto floquet = ComputeFloquetMultipliers( reducedFamily(parameters[i]), data.dataSets[i].points, dt );
		out << parameters[i] << ",";
		for(int64_t j=0;j<floquet.size();++j)
			out << floquet[j].real() << "," << floquet[j].imag() << ",";
		out << endl;
	}
}

void OriginalFloquet( const Options& options ) {

	RosslerFamily rossler;

	auto parameters = ParameterList( 4.0, 8.8, 21 );

	cout << "Generating data..." << endl;
	DataGenerator<RosslerFamily> dataGenerator;
	dataGenerator.initial = Vector3d(7.0,0.0,0.5);
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 100;
	dataGenerator.dtMax = 0.001;
	dataGenerator.print = 10000;
	
	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );

	vector<double> periods(parameters.size());
	for(size_t i=0;i<parameters.size();++i) {
		periods[i] = DetectPeriod(
			data.dataSets[i].points,
			VectorXd(rossler(parameters[i])(data.dataSets[i][0])),
			dataGenerator.tInterval / (dataGenerator.print-1),
			0.001
		);
		cout << periods[i] << endl;
	}

	data = dataGenerator.GenerateDataSystem( parameters, periods, options.numThreads );
	data.WriteDataSetsCSV("output/orig",".csv");
	
	cout << "Computing Floquet multipliers..." << endl;
	ofstream out("output/floquet.csv");
	for(size_t i=0;i<parameters.size();++i) {
		double dt = periods[i] / (dataGenerator.print-1);
		auto floquet = ComputeFloquetMultipliers( rossler(parameters[i]), data.dataSets[i].points, dt );
		out << parameters[i] << ",";
		for(int64_t j=0;j<floquet.size();++j)
			out << floquet[j].real() << "," << floquet[j].imag() << ",";
		out << endl;
	}
}

void ComputeReduced( const Options& options ) {

	RosslerFamily rossler;

	auto parameters = ParameterList( 4.0, 8.8, 21 );

	cout << "Generating data..." << endl;
	DataGenerator<RosslerFamily> dataGenerator;
	dataGenerator.initial = Vector3d(7.0,0.0,0.5);
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 100;
	dataGenerator.dtMax = 0.001;
	dataGenerator.print = 1000;
	
	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );
	data.WriteDataSetsCSV("output/orig",".csv");

	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData( rossler, data, MatrixXd::Identity(3,3), options.numThreads )
	           .WritePointsCSV( "output/p", "-points.csv" )
	           .WriteVectorsCSV( "output/p", "-vectors.csv" )
	           .WriteDerivativesCSV( "output/p", "-derivatives.csv" );

	cout << "Computing Reduced Family..." << endl;
	
	RBFFamilyProducer<RBFType> producer( options.numRBFs );
	producer.boxScale = 1.6;
	auto reducedFamily = producer.BruteForce( reducedData,
											  data.parameters,
											  options.numIterations,
											  options.numThreads );
	
	cout << "Total Cost = " << producer.ComputeTotalCost( reducedFamily, reducedData, data.parameters ) << endl;
	
	reducedFamily.family.WriteCSV("output/reduced.csv");

	cout << "Generating Reduced data..." << endl;
	auto rdataGenerator = MakeDataGenerator( reducedFamily );
	rdataGenerator.MatchSettings( dataGenerator );
	rdataGenerator.tStart = 0.0;

	DataSystem rdata = rdataGenerator.GenerateUsingInitials( parameters, reducedData, options.numThreads );
	rdata.WriteDataSetsCSV("output/rdata",".csv");

	Compare( reducedData, rdata );

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
		BifurcationDiagramGenerator<decltype(reducedFamily)> bifurcationGenerator;
		bifurcationGenerator.tStart = 500.0;
		bifurcationGenerator.tInterval = 500.0;
		bifurcationGenerator.dt = 0.001;
		bifurcationGenerator.dtMax = 0.001;
		bifurcationGenerator.pMin = parameters.front();
		bifurcationGenerator.pMax = parameters.back();
		bifurcationGenerator.pCount = 1280;
		bifurcationGenerator.initial = dataGenerator.initial;
		bifurcationGenerator.Generate( reducedFamily,
			[]( const VectorXd& x, const VectorXd& y ) { return x[1] > 0.0 && y[1] < 0.0; },
			[]( const VectorXd& x, const VectorXd& y ) { return (x[0]+y[0])*0.5; },
			options.numThreads ).WriteBitmap( "output/bifurcation-red.bmp", 720 );
	}
}
/*
void SimulateReduced( const Options& options ) {

	RosslerFamily rossler;

	auto parameters = ParameterList( 4.0, 8.8, 21 );

	cout << "Generating data..." << endl;
	DataGenerator<RosslerFamily> dataGenerator;
	dataGenerator.initial = Vector3d(7.0,0.0,0.5);
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 100;
	dataGenerator.dtMax = 0.001;
	dataGenerator.print = 1000;

	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );
	data.WriteDataSetsCSV("output/orig",".csv");

	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData( rossler, data, MatrixXd::Identity(3,3), options.numThreads );

	RBFFamily<RBFType> reducedFamily;
	reducedFamily.ReadText("reduced.txt");

	RBFFamilyProducer<RBFType> producer( reducedFamily.nRBFs );

	cout << "Total Cost = " << producer.ComputeTotalCost( reducedFamily, reducedData, data.parameters ) << endl;

	producer.Fit( reducedFamily, reducedData, data.parameters );
	
	cout << "Total Cost = " << producer.ComputeTotalCost( reducedFamily, reducedData, data.parameters ) << endl;
	
	reducedFamily.WriteCSV("output/reduced.csv");

	cout << "Generating Reduced data..." << endl;
	DataGenerator<RBFFamily<RBFType>> rdataGenerator( reducedFamily );
	rdataGenerator.MatchSettings( dataGenerator );
	rdataGenerator.tStart = 0.0;

	DataSystem rdata = rdataGenerator.GenerateUsingInitials( parameters, reducedData, options.numThreads );
	rdata.WriteDataSetsCSV("output/rdata",".csv");

	Compare( reducedData, rdata );

	BifurcationDiagramGenerator<RBFFamily<RBFType>> bifurcationGenerator;
	bifurcationGenerator.tStart = 500.0;
	bifurcationGenerator.tInterval = 500.0;
	bifurcationGenerator.dt = 0.001;
	bifurcationGenerator.dtMax = 0.001;
	bifurcationGenerator.pMin = parameters.front();
	bifurcationGenerator.pMax = parameters.back();
	bifurcationGenerator.pCount = 1280;
	bifurcationGenerator.initial = dataGenerator.initial;
	bifurcationGenerator.Generate( reducedFamily,
		[]( const VectorXd& x, const VectorXd& y ) { return x[1] > 0.0 && y[1] < 0.0; },
		[]( const VectorXd& x, const VectorXd& y ) { return (x[0]+y[0])*0.5; },
		options.numThreads ).WriteBitmap( "output/bifurcation-red.bmp", 720 );

}
*/

void Test( const Options& options ) {

	RosslerHighModel::GenerateA(10);

	RosslerHighFamily rosslerHigh(10);

	auto parameters = ParameterList( 4.0, 8.8, 21 );

	cout << "Generating data..." << endl;
	DataGenerator<RosslerHighFamily> dataGenerator(rosslerHigh);
	dataGenerator.initial.setRandom(rosslerHigh.stateDim);
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 100;
	dataGenerator.dtMax = 0.001;
	dataGenerator.print = 1000;
	
	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );

	cout << "Computing secants..." << endl;
	vector<Secants> secants = ComputeSecants( data, 10.0, options.numThreads );

	ProjSecant projSecant( 3 );

	projSecant.ComputeInitial( data )      // Compute initial condition
	          .Find( secants )             // Optimize over Grassmannian
	          .AnalyseSecants( secants )   // Print some statistics
	          .WriteBinary("output/projection.bin")
	          .WriteCSV("output/projection.csv");

	secants = vector<Secants>();

	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData( rosslerHigh, data, projSecant.W, options.numThreads )
	           .WritePointsCSV( "output/p", "-points.csv" )
	           .WriteVectorsCSV( "output/p", "-vectors.csv" )
	           .WriteDerivativesCSV( "output/p", "-derivatives.csv" );

	cout << "Computing Reduced Family..." << endl;
	
	ParameterMapProducer<RosslerFamily> producer;
	RosslerFamily rossler;

	AffineXd parameterMap = producer.SolveSVD( rossler, reducedData, parameters );

	cout << (parameterMap.linear - MatrixXd::Identity(1,1)).norm() << endl;
	cout << (parameterMap.translation - VectorXd::Zero(1)).norm() << endl;

	vector<VectorXd> newParams;
	newParams.reserve( parameters.size() );
	
	transform( begin(parameters), end(parameters), back_inserter(newParams), parameterMap );
	
	cout << "Generating data..." << endl;
	DataGenerator<RosslerFamily> rdataGenerator;
	rdataGenerator.initial = Vector3d(1.0,0.0,0.5);
	rdataGenerator.tStart = 1000;
	rdataGenerator.tInterval = 100;
	rdataGenerator.dtMax = 0.001;
	rdataGenerator.print = 1000;

	DataSystem rdata = rdataGenerator.GenerateDataSystem( newParams, options.numThreads );
	rdata.WriteDataSetsCSV("output/rdata",".csv");

}

void Test2( const Options& options ) {
	RosslerFamily rossler;

	auto parameters = ParameterList( 4.0, 8.8, 21 );

	cout << "Generating data..." << endl;
	DataGenerator<RosslerFamily> dataGenerator;
	dataGenerator.initial = Vector3d(7.0,0.0,0.5);
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 100;
	dataGenerator.dtMax = 0.001;
	dataGenerator.print = 1000;
	
	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );
	data.WriteDataSetsCSV("output/orig",".csv");

	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData( rossler, data, MatrixXd::Identity(3,3), options.numThreads )
	           .WritePointsCSV( "output/p", "-points.csv" )
	           .WriteVectorsCSV( "output/p", "-vectors.csv" )
	           .WriteDerivativesCSV( "output/p", "-derivatives.csv" );

	cout << "Computing Reduced Family..." << endl;
	
	RBFFamilyProducer<RBFType> producer( options.numRBFs );
	producer.boxScale = 1.6;
	auto reducedFamily = producer.BruteForce( reducedData,
											  data.parameters,
											  options.numIterations,
											  options.numThreads );
	
	cout << "Total Cost = " << producer.ComputeTotalCost( reducedFamily, reducedData, data.parameters ) << endl;
}

