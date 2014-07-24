#include <fstream>
#include <DRDSP/projection/proj_pod.h>

using namespace std;
using namespace DRDSP;

ProjPOD& ProjPOD::Find( const DataSet& data ) {

	MatrixXd dataMatrix;
	dataMatrix.setZero(data.dimension,data.points.size());

	uint32_t k = 0;
	for(size_t j=0;j<data.points.size();++j)
		dataMatrix.col(k++) = data.points[j];

	JacobiSVD<MatrixXd> svd(dataMatrix,ComputeThinU);
	singularValues = svd.singularValues();

	W = svd.matrixU().leftCols(targetDimension);
	return *this;
}

ProjPOD& ProjPOD::Find( const DataSystem& data ) {

	size_t totalPoints = 0;
	for(uint32_t i=0;i<data.numParameters;++i)
		totalPoints += data.dataSets[i].points.size();

	MatrixXd dataMatrix;
	dataMatrix.setZero( data.dimension, totalPoints );

	size_t k = 0;
	for(uint32_t i=0;i<data.numParameters;++i)
		for(size_t j=0;j<data.dataSets[i].points.size();++j)
			dataMatrix.col(k++) = data.dataSets[i].points[j];

	JacobiSVD<MatrixXd> svd(dataMatrix,ComputeThinU);
	singularValues = svd.singularValues();

	W = svd.matrixU().leftCols(targetDimension);
	return *this;
}

const ProjPOD& ProjPOD::Write( const char* filename ) const {
	stringstream outfn;
	outfn.str("");
	outfn << filename << "-singularValues.csv";

	ofstream out(outfn.str());
	out.precision(16);
	for(int64_t i=0;i<singularValues.size();++i) {
		out << singularValues[i] << endl;
	}
	out.close();

	outfn.str("");
	outfn << filename << "-W.csv";
	out.open(outfn.str());
	out.precision(16);
	for(int64_t i=0;i<W.rows();++i) {
		for(int64_t j=0;j<W.cols();++j)
			out << W(i,j) << ",";
		out << endl;
	}
	out.close();

	outfn.str("");
	outfn << filename << "-W.bin";
	out.open(outfn.str(),ios::binary);
	for(int64_t i=0;i<W.rows();++i)
		for(int64_t j=0;j<W.cols();++j)
			out.write((char*)&W(i,j),sizeof(double));
	out.close();
	return *this;
}
