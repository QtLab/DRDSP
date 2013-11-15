#include <fstream>
#include <DRDSP/dynamics/model_rbf.h>

using namespace std;
using namespace DRDSP;

ModelRBF::ModelRBF() : weights(nullptr), rbfs(nullptr), numRBFs(0) {}

ModelRBF::ModelRBF( uint16_t dim, uint16_t nRBFs ) : weights(nullptr), rbfs(nullptr), numRBFs(0) {
	Create(dim,nRBFs);
}

ModelRBF::ModelRBF( const ModelRBF& rhs ) {
	Create(rhs.dimension,rhs.numRBFs);
	linear = rhs.linear;
	for(uint16_t i=0;i<numRBFs;i++) {
		weights[i] = rhs.weights[i];
		rbfs[i] = rhs.rbfs[i];
	}
}

ModelRBF::ModelRBF( ModelRBF&& rhs ) {
	dimension = rhs.dimension;
	numRBFs = rhs.numRBFs;
	linear = rhs.linear;
	weights = rhs.weights;
	rbfs = rhs.rbfs;
	rhs.numRBFs = 0;
	rhs.weights = nullptr;
	rhs.rbfs = nullptr;
}

ModelRBF& ModelRBF::operator=( const ModelRBF& rhs ) {
	Destroy();
	Create(rhs.dimension,rhs.numRBFs);
	linear = rhs.linear;
	for(uint16_t i=0;i<numRBFs;i++) {
		weights[i] = rhs.weights[i];
		rbfs[i] = rhs.rbfs[i];
	}
	return *this;
}

ModelRBF& ModelRBF::operator=( ModelRBF&& rhs ) {
	if( this != &rhs ) {
		Destroy();
		dimension = rhs.dimension;
		numRBFs = rhs.numRBFs;
		linear = rhs.linear;
		weights = rhs.weights;
		rbfs = rhs.rbfs;
		rhs.numRBFs = 0;
		rhs.weights = nullptr;
		rhs.rbfs = nullptr;
	}
	return *this;
}

void ModelRBF::Create( uint16_t dim, uint16_t nRBFs ) {
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

void ModelRBF::SetRBFType( const Function& f ) {
	for(uint16_t i=0;i<numRBFs;i++) {
		rbfs[i].function = &f;
	}
}

VectorXd ModelRBF::VectorField( const VectorXd& x ) {
	VectorXd sum = linear * x;
	for(uint16_t i=0;i<numRBFs;i++)
		sum += weights[i] * rbfs[i](x);
	return sum;
}

MatrixXd ModelRBF::Partials( const VectorXd &x ) {
	MatrixXd sum = linear;
	for(uint16_t i=0;i<numRBFs;i++)
		sum += weights[i] * rbfs[i].Derivative(x).transpose();
	return sum;
}

void ModelRBF::SetCentresRandom( const AABB& box ) {
	double rnd;
	VectorXd diff = box.bMax - box.bMin;
	for(uint16_t i=0;i<numRBFs;i++)
		for(uint16_t j=0;j<dimension;j++) {
			rnd = (double)rand()/RAND_MAX;
			rbfs[i].centre(j) = box.bMin(j) + diff(j) * rnd;
		}
}

void ModelRBF::LoadCentresText( const char* filename ) {
	ifstream in(filename);
	if( !in ) return;

	for(uint16_t k=0;k<numRBFs;k++)
		for(uint16_t j=0;j<dimension;j++)
			in >> rbfs[k].centre(j);
	in.close();
}

void ModelRBF::LoadCentresBinary( const char* filename ) {
	ifstream in(filename);
	if( !in ) return;

	for(uint16_t k=0;k<numRBFs;k++)
		in.read( (char*)&rbfs[k].centre(0), sizeof(double)*dimension );
	in.close();
}

void ModelRBF::WriteCSV( const char *filename ) const {
	
	ofstream out;
	out.open(filename);
	out.precision(16);
	out << dimension << "," << numRBFs << endl;
	for(uint16_t i=0;i<dimension;i++) {	
		for(uint16_t j=0;j<dimension;j++)
			out << linear(i,j) << ",";
		out << endl;
	}
	out << endl;
	for(uint16_t i=0;i<dimension;i++) {	
		for(uint16_t j=0;j<numRBFs;j++)
			out << weights[j](i) << ",";
		out << endl;
	}
	out << endl;
	for(uint16_t i=0;i<dimension;i++) {	
		for(uint16_t j=0;j<numRBFs;j++)
			out << rbfs[j].centre(i) << ",";
		out << endl;
	}
	out.close();
	
}

