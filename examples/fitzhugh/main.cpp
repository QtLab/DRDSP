#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/rbf_family_producer.h>
#include <DRDSP/dynamics/data_generator.h>
#include <DRDSP/misc.h>
#include <DRDSP/projection/inverse.h>
#include "fitzhugh_nagumo.h"

using namespace std;
using namespace DRDSP;

struct Options {
	uint32_t targetDimension, numRBFs, numIterations, numThreads;

	Options() : targetDimension(2), numRBFs(30), numIterations(1000), numThreads(3) {}
	
	Options( int argc, char** argv ) : Options() {
		if( argc >= 2 ) targetDimension = (uint32_t)atoi(argv[1]);
		if( argc >= 3 )         numRBFs = (uint32_t)atoi(argv[2]);
		if( argc >= 4 )   numIterations = (uint32_t)atoi(argv[3]);
		if( argc >= 5 )      numThreads = (uint32_t)atoi(argv[4]);
	}
};

typedef EquiRBFZ2<ThinPlateSpline> RBFType;

void ComputeReduced( const Options& options );

int main( int argc, char** argv ) {
	Options options(argc,argv);

	ComputeReduced( options );

	cout << "Press enter to continue . . . "; cin.get();
}

void ComputeReduced( const Options& options ) {

	FitzHughNagumoFamily fitzHughNagumoFamily;

	auto parameters = ParameterList( 0.5, 4.0, 12 );

	cout << "Generating data..." << endl;
	DataGenerator<FitzHughNagumoFamily> dataGenerator;
	dataGenerator.initial.setRandom( fitzHughNagumoFamily.stateDim );
	dataGenerator.initial.array() -= 0.5;
	dataGenerator.initial *= 2.0;
	dataGenerator.tStart = 100.0;
	dataGenerator.tInterval = 20.0;
	dataGenerator.print = 1000;
	dataGenerator.dtMax = 0.001;

	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );

	cout << "Computing secants..." << endl;
	vector<Secants> secants = ComputeSecants( data, 10.0, options.numThreads );

	ProjSecant projSecant( options.targetDimension );

	projSecant.ComputeInitial( data )      // Compute initial condition
	          .Find( secants )             // Optimize over Grassmannian
	          .AnalyseSecants( secants )   // Print some statistics
	          .WriteBinary("output/projection.bin")
	          .WriteCSV("output/projection.csv");

	secants = vector<Secants>();

	auto inverse = ComputeInverse( projSecant.W, data );
	cout << ComputeInverseCost( inverse, data ) << endl;

	cout << endl << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData( fitzHughNagumoFamily, data, projSecant.W, options.numThreads )
	           .WritePointsCSV("output/p","-points.csv")
	           .WriteVectorsCSV("output/p","-vectors.csv");

	cout << endl << "Computing Reduced Family..." << endl;
	RBFFamilyProducer<RBFType> producer( options.numRBFs );
	auto reducedFamily = producer.BruteForce( reducedData,
											  data.parameters,
											  options.numIterations,
											  options.numThreads );

	cout << "Total Cost = "
	     << producer.ComputeTotalCost( reducedFamily, reducedData, data.parameters )
	     << endl;

	reducedFamily.family.WriteCSV("output/reduced.csv");

	cout << "Simulating the reduced family..." << endl;
	auto rdataGenerator = MakeDataGenerator( reducedFamily );
	rdataGenerator.MatchSettings( dataGenerator );
	rdataGenerator.tStart = 0.0;
	DataSystem rdata = rdataGenerator.GenerateUsingInitials( parameters, reducedData, options.numThreads );
	rdata.WriteCSV("output/rdata",".csv");

	Compare( reducedData, rdata );
	ComparePeriods( reducedData, rdata, dataGenerator.tInterval / (dataGenerator.print-1), 1.0 );
}
