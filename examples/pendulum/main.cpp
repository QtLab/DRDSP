#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/rbf_family_producer.h>
#include <DRDSP/dynamics/data_generator.h>
#include <DRDSP/misc.h>
#include "pendulum.h"
#include <DRDSP/dynamics/bifurcation.h>

using namespace std;
using namespace DRDSP;

struct Options {
	uint32_t targetDimension = 3,
	         numRBFs = 100,
	         numIterations = 120,
	         numThreads = 4;
	
	Options( int argc, char** argv ) {
		if( argc >= 2 ) targetDimension = (uint32_t)atoi(argv[1]);
		if( argc >= 3 )         numRBFs = (uint32_t)atoi(argv[2]);
		if( argc >= 4 )   numIterations = (uint32_t)atoi(argv[3]);
		if( argc >= 5 )      numThreads = (uint32_t)atoi(argv[4]);
	}
};

using RBFType = RBF<PolyharmonicSpline<4>>;

int main( int argc, char** argv ) {

	Options options(argc,argv);

	FamilyEmbedded<PendulumFamily,FlatEmbedding> pendulum;

	auto parameters = ParameterList( 1.886, 1.894, 8 );

	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator<PendulumFamily,PendulumSolver> dataGenerator;
	dataGenerator.initial.setZero( pendulum.family.stateDim );
	dataGenerator.initial[0] = -0.81;
	dataGenerator.initial[1] = 0.14;
	dataGenerator.initial[3] = 0.74;
	dataGenerator.tStart = 1000;
	dataGenerator.tInterval = 1000;
	dataGenerator.print = 2000;
	dataGenerator.dtMax = 0.001;

	DataSystem data = dataGenerator.GenerateDataSystem( parameters, options.numThreads );
	data.WriteCSV("output/orig",".csv");
	
	cout << "Embedding data..." << endl;
	DataSystem dataEmbedded = EmbedData( pendulum.embedding, data, options.numThreads );

	cout << "Computing secants..." << endl;
	vector<Secants> secants = ComputeSecants( dataEmbedded, 10.0, options.numThreads );

	cout << "Finding projection..." << endl;
	ProjSecant projSecant( options.targetDimension );
	
	projSecant.ComputeInitial( dataEmbedded )   // Compute initial condition
	          .Find( secants )                  // Optimize over Grassmannian
	          .AnalyseSecants( secants )        // Print some statistics
	          .WriteBinary("output/projection.bin")
	          .WriteCSV("output/projection.csv");

	secants = vector<Secants>();

/*
	BifurcationDiagramGenerator<PendulumFamily,PendulumSolver> bifurcationGenerator;
	bifurcationGenerator.tStart = 1000.0;
	bifurcationGenerator.tInterval = 1000.0;
	bifurcationGenerator.dt = 0.001;
	bifurcationGenerator.dtMax = 0.001;
	bifurcationGenerator.pMin = parameters.front();
	bifurcationGenerator.pMax = parameters.back();
	bifurcationGenerator.pCount = 1280;
	bifurcationGenerator.initial = dataGenerator.initial;
	bifurcationGenerator.Generate( pendulum.family,
		[&]( const VectorXd& x, const VectorXd& y ) {
			VectorXd px = projSecant.W.adjoint() * pendulum.embedding(x);
			VectorXd py = projSecant.W.adjoint() * pendulum.embedding(y);
			return px[2] < 0.0 && py[2] > 0.0;
		},
		[&]( const VectorXd& x, const VectorXd& y ) {
			VectorXd px = projSecant.W.adjoint() * pendulum.embedding(x);
			VectorXd py = projSecant.W.adjoint() * pendulum.embedding(y);
			return (px[1]+py[1])*0.5;
		},
		options.numThreads ).WriteBitmap( "output/bifurcation-orig.bmp", 720 );
*/


	cout << "Computing Reduced Data..." << endl;
	ReducedDataSystem reducedData;
	reducedData.ComputeDataEmbedded( pendulum, data, projSecant.W, options.numThreads )
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

	cout << "Simulating the reduced model..." << endl;
	auto rdataGenerator = MakeDataGenerator( reducedFamily );
	rdataGenerator.MatchSettings(dataGenerator);
	rdataGenerator.tStart = 0.0;
	DataSystem rdata = rdataGenerator.GenerateUsingInitials( parameters, reducedData, options.numThreads );
	rdata.WriteCSV("output/rdata",".csv");

	Compare( reducedData, rdata );

	cout << "Press enter to continue . . . "; cin.get();
}
