#include <iostream>
#include <fstream>
#include <vector>
#include <DRDSP/misc.h>
#include <DRDSP/projection/proj_secant.h>
#include <DRDSP/optimization/gradient_descent.h>
#include <DRDSP/optimization/conjugate_gradient.h>
#include <DRDSP/geometry/grassmannian.h>

using namespace std;
using namespace DRDSP;


double SecantCostFunction::operator()( const MatrixXd& X ) const {
	
	if( secants.count == 0 ) return 1.0;

	double sum = 0.0;

	if( secants.weights.size() > 0 ) {
		uint32_t sumWeights = accumulate( secants.weights, uint32_t(0) );
		for(size_t j=0;j<secants.count;++j) {
			sum += (double)secants.weights[j] / ( X.adjoint() * secants.GetSecant(j) ).norm();
		}
		sum *= 1.0 / sumWeights;
	} else {
		for( const auto& k : secants.secants ) {
			sum += 1.0 / ( X.adjoint() * k ).norm();
		}
		sum *= 1.0 / secants.count;
	}
	return sum;
}

double SecantCostFunctionMulti::operator()( const MatrixXd& X ) const {
	
	double sum = 0.0;

	for(size_t i=0;i<N;++i) {
		sum += SecantCostFunction(secants[i])(X);
	}

	return sum / N;
}

MatrixXd SecantCostGradient::operator()( const MatrixXd& X ) const {

	MatrixXd sum;
	sum.setZero(X.rows(),X.cols());

	if( secants.count == 0 ) return sum;

	double projectedLength;
	VectorXd secant, projectedSecant;

	if( secants.weights.size() > 0 ) {

		uint32_t sumWeights = accumulate( secants.weights, uint32_t(0) );

		for(size_t j=0;j<secants.count;++j) {
			secant = secants.GetSecant(j);
			projectedSecant = X.adjoint() * secant;
			projectedLength = projectedSecant.norm();
			sum += (((double)secants.weights[j] / ( projectedLength * projectedLength * projectedLength )) * secant) * projectedSecant.transpose();
		}

		sum *= 1.0 / sumWeights;

	} else {
		for(size_t j=0;j<secants.count;++j) {
			secant = secants.GetSecant(j);
			projectedSecant = X.adjoint() * secant;
			projectedLength = projectedSecant.norm();
			sum += ((1.0 / ( projectedLength * projectedLength * projectedLength )) * secant) * projectedSecant.transpose();
		}
		sum *= 1.0 / secants.count;
	}

	return Grassmannian::HorizontalComponent( X, -sum );
}

MatrixXd SecantCostGradientMulti::operator()( const MatrixXd& X ) const {
	
	MatrixXd sum = SecantCostGradient(secants[0])(X);

	for(size_t i=1;i<N;++i)
		sum += SecantCostGradient(secants[i])(X);

	return sum / (double)N;
}

ProjSecant::ProjSecant() : ProjSecant(2) {}

ProjSecant::ProjSecant( uint32_t targetDimension ) :
	targetMinProjectedLength(0.5),
	targetDimension(targetDimension),
	maxIterations(100)
{}

ProjSecant& ProjSecant::Find( const Secants& secants ) {
	SecantCostFunction S(secants);
	SecantCostGradient gradS(secants);

	ConjugateGradient<Grassmannian::Geodesic,SecantCostFunction,SecantCostGradient> optimiziation( S, gradS );

	optimiziation.maxSteps = maxIterations;
	optimiziation.lineSearch.alpha = 2.0;
	optimiziation.Optimize( W );
	return *this;
}

ProjSecant& ProjSecant::Find( const Secants* secants, size_t N ) {
	SecantCostFunctionMulti S(secants,N);
	SecantCostGradientMulti gradS(secants,N);

	ConjugateGradient<Grassmannian::Geodesic,SecantCostFunctionMulti,SecantCostGradientMulti> optimiziation( S, gradS );

	optimiziation.maxSteps = maxIterations;
	optimiziation.lineSearch.alpha = 2.0;
	optimiziation.Optimize( W );
	return *this;
}

ProjSecant& ProjSecant::Find( const vector<Secants>& secants ) {
	return Find( secants.data(), secants.size() );
}

ProjSecant& ProjSecant::ComputeInitial( const DataSet& data ) {
	uint32_t n = data.dimension;

	vector<double> maxVal(n);
	vector<double> minVal(n);
	vector<double> spread(n);
	double val, bigVal;
	uint32_t bigAxis;

	for(uint32_t k=0;k<n;++k)
		maxVal[k] = minVal[k] = data.points[0](0);

	for(size_t j=0;j<data.points.size();++j)
		for(uint32_t k=0;k<n;++k) {
			val = data.points[j](k);
			if( val > maxVal[k] )
				maxVal[k] = val;
			if( val < minVal[k] )
				minVal[k] = val;
		}

	for(uint32_t k=0;k<n;++k)
		spread[k] = maxVal[k] - minVal[k];


	W.setZero(n,targetDimension);

	cout << "Initial Condition: ( ";

	for(uint32_t i=0;i<targetDimension;++i) {
		bigVal = 0.0;
		bigAxis = 0;
		for(uint32_t k=0;k<n;++k) {
			if( spread[k] > bigVal ) {
				bigVal = spread[k];
				bigAxis = k;
			}
		}
		spread[bigAxis] = 0.0;
		cout << bigAxis;
		if( i < targetDimension - 1 )  cout << ", ";
		W(bigAxis,i) = 1.0;
	}
	cout << " )" << endl;
	return *this;
}

ProjSecant& ProjSecant::ComputeInitial( const DataSystem& data ) {
	uint32_t n = data.dimension;

	vector<double> maxVal(n);
	vector<double> minVal(n);
	vector<double> spread(n);
	double val, bigVal;
	uint32_t bigAxis;

	for(uint32_t k=0;k<n;++k)
		maxVal[k] = minVal[k] = data.dataSets[0].points[0](k);

	for(uint32_t i=0;i<data.numParameters;++i)
		for(size_t j=0;j<data.dataSets[i].points.size();++j)
			for(uint32_t k=0;k<n;++k) {
				val = data.dataSets[i].points[j](k);
				if( val > maxVal[k] )
					maxVal[k] = val;
				if( val < minVal[k] )
					minVal[k] = val;
			}

	for(uint32_t k=0;k<n;++k)
		spread[k] = maxVal[k] - minVal[k];

	W.setZero(n,targetDimension);

	cout << "Initial Condition: ( ";

	for(uint32_t i=0;i<targetDimension;++i) {
		bigVal = 0.0;
		bigAxis = 0;
		for(uint32_t k=0;k<n;++k) {
			if( spread[k] > bigVal ) {
				bigVal = spread[k];
				bigAxis = k;
			}
		}
		spread[bigAxis] = 0.0;
		cout << bigAxis;
		if( i < targetDimension - 1 )  cout << ", ";
		W(bigAxis,i) = 1.0;
	}
	cout << " )" << endl;
	return *this;
}

const ProjSecant& ProjSecant::AnalyseSecants( const Secants& secants ) const {
	return AnalyseSecants(&secants,1);
}

const ProjSecant& ProjSecant::AnalyseSecants( const Secants* secants, size_t N ) const {
	double xMin = 1.0, xMax = 0.0, xMean = 0.0, total = 0.0, len;
	size_t numSecants = 0;
	for(size_t i=0;i<N;++i) {
		for(size_t j=0;j<secants[i].count;++j) {
			len = ( W.adjoint() * secants[i].GetSecant(j) ).norm();
			if( len < xMin ) xMin = len;
			if( len > xMax ) xMax = len;
			total += len;
		}
		numSecants += secants[i].count;
	}
	xMean = total / numSecants;
	cout << endl << "Projected Lengths: Range = [ " << xMin << ", " << xMax << " ], Mean = " << xMean << endl;
	return *this;
}

const ProjSecant& ProjSecant::AnalyseSecants( const vector<Secants>& secants ) const {
	return AnalyseSecants(secants.data(),secants.size());
}

const ProjSecant& ProjSecant::WriteCSV( const char* filename ) const {
	ofstream out(filename);
	out.precision(16);
	for(int64_t i=0;i<W.rows();++i) {
		for(int64_t j=0;j<W.cols();++j)
			out << W(i,j) << ",";
		out << endl;
	}
	return *this;
}

const ProjSecant& ProjSecant::WriteBinary( const char* filename ) const {
	ofstream out(filename,ios::binary);
	for(int64_t i=0;i<W.rows();++i)
		for(int64_t j=0;j<W.cols();++j)
			out.write((char*)&W(i,j),sizeof(double));
	return *this;
}

bool ProjSecant::ReadBinary( const char* filename ) {
	ifstream in(filename,ios::binary);
	if( !in ) {
		cout << "ProjSecant::Read : file error" << endl;
		return false;
	}

	in.seekg(0, ios::end);
	if( (size_t)in.tellg() < sizeof(double)*W.size() ) {
		cout << "ProjSecant::Read : insufficient data" << endl;
		return false;
	}
	in.seekg(0, ios::beg);

	for(int64_t i=0;i<W.rows();++i)
		for(int64_t j=0;j<W.cols();++j)
			in.read((char*)&W(i,j),sizeof(double));

	return true;
}
