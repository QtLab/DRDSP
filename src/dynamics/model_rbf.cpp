#include <fstream>
#include <DRDSP/dynamics/model_rbf.h>

using namespace std;
using namespace DRDSP;

ModelRBF::ModelRBF() : weights(nullptr), rbfs(nullptr), dimension(0), numRBFs(0) {}

ModelRBF::ModelRBF( uint32_t dim, uint32_t nRBFs ) : weights(nullptr), rbfs(nullptr), dimension(0), numRBFs(0) {
	Create(dim,nRBFs);
}

void ModelRBF::Create( uint32_t dim, uint32_t nRBFs ) {
	dimension = dim;
	numRBFs = nRBFs;
	weights = new VectorXd [numRBFs];
	rbfs = new RadialFunction [numRBFs];
	linear.setZero(dimension,dimension);
	for(uint16_t i=0;i<numRBFs;i++) {
		weights[i].setZero(dimension);
		rbfs[i].centre.setZero(dimension);
	}
}

void ModelRBF::SetRBFType( const Function& f ) {
	for(uint16_t i=0;i<numRBFs;i++) {
		rbfs[i].function = &f;
	}
}

void ModelRBF::Destroy() {
	if(weights) {
		delete[] weights;
		weights = nullptr;
	}
	if(rbfs) {
		delete[] rbfs;
		rbfs = nullptr;
	}
	numRBFs = 0;
}

VectorXd ModelRBF::VectorField( const VectorXd& x ) const {
	VectorXd sum = linear * x;
	for(uint16_t i=0;i<numRBFs;i++)
		sum += weights[i] * rbfs[i](x);
	return sum;
}

MatrixXd ModelRBF::VectorFieldDerivative( const VectorXd &x ) const {
	MatrixXd sum = linear;
	for(uint32_t i=0;i<numRBFs;i++)
		sum += weights[i] * rbfs[i].Derivative(x).transpose();
	return sum;
}

void ModelRBF::SetCentresRandom( const VectorXd& minBounds, const VectorXd& maxBounds ) {
	double rnd;
	VectorXd diff = maxBounds - minBounds;
	for(uint32_t i=0;i<numRBFs;i++)
		for(uint32_t j=0;j<dimension;j++) {
			rnd = (double)rand()/RAND_MAX;
			rbfs[i].centre(j) = minBounds(j) + diff(j) * rnd;
		}
}

void ModelRBF::LoadCentresText( const char* filename ) {
	ifstream in(filename);
	if( !in ) return;

	for(uint k=0;k<numRBFs;k++)
		for(uint j=0;j<dimension;j++)
			in >> rbfs[k].centre(j);
	in.close();
}

void ModelRBF::LoadCentresBinary( const char* filename ) {
	ifstream in(filename);
	if( !in ) return;

	for(uint32_t k=0;k<numRBFs;k++)
		in.read( (char*)&rbfs[k].centre(0), sizeof(double)*dimension );
	in.close();
}


