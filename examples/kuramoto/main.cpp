#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/rbf_family_producer.h>
#include <DRDSP/dynamics/data_generator.h>
#include "kuramoto.h"

using namespace std;
using namespace DRDSP;

struct Options {
	uint32_t numIterations, maxPoints, targetDimension, numRBFs;

	Options() : numIterations(500), maxPoints(0), targetDimension(4), numRBFs(20) {}
	
	Options( int argc, char** argv ) : Options() {
		if( argc >= 2 ) targetDimension = (uint32_t)atoi(argv[1]);
		if( argc >= 3 ) numRBFs = (uint32_t)atoi(argv[2]);
		if( argc >= 4 ) maxPoints = (uint32_t)atoi(argv[3]);
		if( argc >= 5 ) numIterations = (uint32_t)atoi(argv[4]);
	}
};

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

	Options options(argc,argv);;

	// The kuramoto example
	FamilyEmbedded<KuramotoAFamily,FlatEmbedding> kuramoto(KuramotoAFamily(100),FlatEmbedding(101));
	
	auto parameters = ParameterList( 1.0, 1.1, 6 );

	// Generate the data
	cout << "Generating data..." << endl;
	DataGenerator<KuramotoAFamily,KuramotoASolver> dataGenerator(kuramoto.family);
	dataGenerator.initial.setRandom();
	dataGenerator.tStart = 50;
	dataGenerator.tInterval = 7;
	dataGenerator.print = 100;
	dataGenerator.dtMax = 0.001;

	DataSystem data = dataGenerator.GenerateDataSystem( parameters, 4 );
			
	// Embed the data
	cout << "Embedding data..." << endl;
	DataSystem dataEmbedded = EmbedData( kuramoto.embedding, data );

	// Pre-compute secants
	cout << "Computing secants..." << endl;
	vector<SecantsPreComputed> secants( dataEmbedded.numParameters );

	for(uint32_t i=0;i<dataEmbedded.numParameters;++i)
		secants[i].ComputeFromData( dataEmbedded.dataSets[i] );

	// Secant culling
	cout << "Culling secants..." << endl;
	vector<SecantsPreComputed> newSecants( dataEmbedded.numParameters );
	for(uint32_t i=0;i<dataEmbedded.numParameters;++i)
		newSecants[i] = secants[i].CullSecantsDegrees( 10.0 );

	secants = vector<SecantsPreComputed>();

	// Find a projection
	cout << "Finding projection..." << endl;
	ProjSecant projSecant;
	projSecant.targetDimension = options.targetDimension;
	projSecant.targetMinProjectedLength = 0.7;

	// Compute initial condition
	// For this particular example, we use a custom initial condition
	projSecant.W.setZero(kuramoto.embedding.eDim,4);
	
	for(uint32_t i=0;i<(kuramoto.embedding.eDim-2)/2;i++) {
		projSecant.W(2*i,0) = 1;
		projSecant.W(2*i+1,1) = 1;
	}
	projSecant.W(kuramoto.embedding.eDim-2,2) = 1;
	projSecant.W(kuramoto.embedding.eDim-1,3) = 1;
	
	projSecant.W.col(0).normalize();
	projSecant.W.col(1).normalize();

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
	reducedData.ComputeDataEmbedded( kuramoto, data, projSecant.W, 4 );
	reducedData.WritePointsCSV("output/p","-points.csv");
	reducedData.WriteVectorsCSV("output/p","-points.csv");

	// Obtain the reduced model
	cout << "Computing Reduced Model..." << endl;
	
	RBFFamilyProducer<RadialType> producer(options.numRBFs);
	auto reducedModel = producer.BruteForce(reducedData,data.parameterDimension,data.parameters,options.numIterations);
	
	cout << "Total Cost = " << producer.ComputeTotalCost(reducedModel,reducedData,data.parameters) << endl;
	
	reducedModel.WriteCSV("output/reduced.csv");
	
	//ModelRBFProducer rbfProducer(options.numRBFs);
	//ModelRBF rbfModel = rbfProducer.BruteForce(reducedData.reducedData[0],options.numIterations);

	//cout << "Total Cost = " << rbfProducer.ComputeTotalCost(rbfModel,reducedData.reducedData[0]) << endl;	
	
	// Generate the data
	cout << "Generating Reduced data..." << endl;
	DataGenerator<RBFFamily<RadialType>> rdataGenerator(reducedModel);
	rdataGenerator.MatchSettings(dataGenerator);
	rdataGenerator.tStart = 0.0;

	DataSystem rdata = rdataGenerator.GenerateUsingInitials( parameters, reducedData, 4 );
	rdata.WriteDataSetsCSV("output/rdata",".csv");

	Compare( reducedData, rdata );

	system("PAUSE");
}
