#include <DRDSP/data/data_set.h>
#include <DRDSP/data/secants.h>
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/dynamics/model_rbf_producer.h>

//#include "pendulum.h"

using namespace DRDSP;

int main() {

	// Create a new data set
	DataSet data;

	// Load data from file
	data.LoadSetBinary("data/p1.8.bin");

	// Pre-compute secants if less than 128 MB
	Secants secants;
	secants.ComputeFromData( data, 1 << 27 );

	// Secant culling
	Secants newSecants = secants.CullSecantsDegrees( 10.0 );

	// Find a projection
	ProjSecant projSecant;
	projSecant.targetDimension = 2;
	projSecant.targetMinProjectedLength = 0.7;

	// Compute initial condition
	projSecant.GetInitial(data);

	// Optimize over Grassmannian
	projSecant.Find(newSecants);

	// Print some statistics
	projSecant.AnalyseSecants(newSecants);

	// Dynamics
	//Pendulum brusselator;
	VectorXd parameter(1);
	parameter(0) = 2.0;

	// Compute projected data
	ReducedData reducedData;
	//reducedData.ComputeData(brusselator,data,parameter,projSecant.W);

	// Obtain the reduced model
	ModelRBFProducer producer;
	//producer.numRBFs = numRBFs;

	// Find a reduced model for the projected data using 100 radial basis functions
	ModelRBF reducedModel = producer.BruteForce(reducedData,100);


	return 0;
}
