#include <fstream>
#include <DRDSP/projection/proj_pod.h>

using namespace std;
using namespace DRDSP;

ProjPOD::ProjPOD() : targetDimension(2) {}

ProjPOD::ProjPOD( uint32_t targetDimension ) : targetDimension(targetDimension) {}

void ProjPOD::Find( const DataSet& data ) {

	dataMatrix.setZero(data.dimension,data.points.size());

	uint32_t k = 0;
	for(uint32_t j=0;j<data.points.size();j++)
		dataMatrix.col(k++) = data.points[j];

	svd.compute(dataMatrix,ComputeThinU);
	
	W = svd.matrixU().block(0,0,data.dimension,targetDimension);
}

void ProjPOD::Find( const DataSystem& data ) {

	size_t totalPoints = 0;
	for(uint16_t i=0;i<data.numParameters;i++)
		totalPoints += data.dataSets[i].points.size();

	dataMatrix.setZero( data.dimension, totalPoints );

	size_t k = 0;
	for(uint16_t i=0;i<data.numParameters;i++)
		for(size_t j=0;j<data.dataSets[i].points.size();j++)
			dataMatrix.col(k++) = data.dataSets[i].points[j];

	svd.compute(dataMatrix,ComputeThinU);
	
	W = svd.matrixU().block(0,0,data.dimension,targetDimension);
}

void ProjPOD::Write( const char* filename ) const {
	stringstream outfn;
	outfn.str("");
	outfn << filename << "-singularValues.csv";

	ofstream out(outfn.str());
	out.precision(16);
	for(int i=0;i<svd.nonzeroSingularValues();i++) {
		out << svd.singularValues()[i] << endl;
	}
	out.close();

	outfn.str("");
	outfn << filename << "-W.csv";
	out.open(outfn.str());
	out.precision(16);
	for(int i=0;i<W.rows();i++) {
		for(int j=0;j<W.cols();j++)
			out << W(i,j) << ",";
		out << endl;
	}
	out.close();

	outfn.str("");
	outfn << filename << "-W.bin";
	out.open(outfn.str(),ios::binary);
	for(int i=0;i<W.rows();i++)
		for(int j=0;j<W.cols();j++)
			out.write((char*)&W(i,j),sizeof(double));
	out.close();
}
