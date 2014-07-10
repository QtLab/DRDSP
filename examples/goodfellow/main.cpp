#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/rbf_family_producer.h>
#include <DRDSP/dynamics/data_generator.h>
#include <DRDSP/misc.h>
#include "goodfellow.h"

using namespace std;
using namespace DRDSP;

struct Options {
	uint32_t targetDimension, numRBFs, numIterations, numThreads;

	Options() : targetDimension(2), numRBFs(100), numIterations(300), numThreads(3) {}
	
	Options( int argc, char** argv ) : Options() {
		if( argc >= 2 ) targetDimension = (uint32_t)atoi(argv[1]);
		if( argc >= 3 )         numRBFs = (uint32_t)atoi(argv[2]);
		if( argc >= 4 )   numIterations = (uint32_t)atoi(argv[3]);
		if( argc >= 5 )      numThreads = (uint32_t)atoi(argv[4]);
	}
};

typedef EquiRBFZ2<Multiquadratic> RBFType; //

int main( int argc, char** argv ) {

	Options options(argc,argv);

	GoodfellowFamily goodfellow(10);

	auto parameters = ParameterList( 4.1, 4.2, 6 );
	
	cout << "Generating data..." << endl;
	DataGenerator<GoodfellowFamily> dataGenerator( goodfellow );
	dataGenerator.initial.setRandom( goodfellow.stateDim );
	dataGenerator.tStart = 100;
	dataGenerator.tInterval = 2;
	dataGenerator.print = 200;
	dataGenerator.dtMax = 0.00001;

	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );
	data.WriteDataSetsCSV("output/orig",".csv");

	cout << "Computing secants..." << endl;
	vector<Secants> secants = ComputeSecants( data, 10.0, options.numThreads );

	cout << "Finding projection..." << endl;
	ProjSecant projSecant( options.targetDimension );

	projSecant.ComputeInitial( data )      // Compute initial condition
	          .Find( secants )             // Optimize over Grassmannian
	          .AnalyseSecants( secants )   // Print some statistics
	          .WriteBinary("output/projection.bin")
	          .WriteCSV("output/projection.csv");

	secants = vector<Secants>();

	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeData( goodfellow, data, projSecant.W, options.numThreads )
	           .WritePointsCSV("output/p","-points.csv")
	           .WriteVectorsCSV("output/p","-vectors.csv");

	cout << "Computing Reduced Family..." << endl;
	
	RBFFamilyProducer<RBFType> producer( options.numRBFs );
	auto reducedFamily = producer.BruteForce( reducedData,
											  data.parameters,
											  options.numIterations,
											  options.numThreads );
	
	cout << "Total Cost = "
	     << producer.ComputeTotalCost( reducedFamily, reducedData, data.parameters )
		 << endl;
	
	reducedFamily.family.WriteCSV("output/reduced.csv");

	cout << "Generating Reduced data..." << endl;
	auto rdataGenerator = MakeDataGenerator( reducedFamily );
	rdataGenerator.MatchSettings( dataGenerator );
	rdataGenerator.tStart = 0.0;

	DataSystem rdata = rdataGenerator.GenerateUsingInitials( parameters, reducedData, options.numThreads );
	rdata.WriteDataSetsCSV("output/rdata",".csv");

	Compare( reducedData, rdata );

	cout << "Press enter to continue . . . "; cin.get();
}

