#include <DRDSP/projection/proj_system.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <DRDSP/optimization/gradient_descent.h>
#include <DRDSP/optimization/conjugate_gradient.h>
#include <DRDSP/geometry/grassmannian.h>
#include <DRDSP/geometry/metric.h>

using namespace std;
using namespace DRDSP;

ProjSecantSystem::ProjSecantSystem( System *sys ) : data(sys), minPLen(1.0), d(3), maxIts(50) {
	n = data->dim;
	Proj = new ProjSingle* [data->numFiles];
	for(uint32_t i=0;i<data->numFiles;i++)
		Proj[i] = new ProjSingle( &data->file[i] );
}

ProjSecantSystem::~ProjSecantSystem() {
	for(uint32_t i=0;i<data->numFiles;i++)
		delete Proj[i];
	delete[] Proj;
}

void ProjSecantSystem::GetInitial() {

	double *maxVal = new double [n];
	double *minVal = new double [n];
	double *spread = new double [n];
	double val, bigVal;
	uint bigAxis;

	for(uint32_t k=0;k<n;k++) {
		maxVal[k] = 0.0;
		minVal[k] = 0.0;
	}

	for( uint32_t i=0;i<data->numFiles;i++ )
		for(uint32_t j=0;j<data->file[i].Points->count;j++)
			for(uint32_t k=0;k<n;k++) {
				val = data->file[i].Points->p[j]->pos(k);
				if( i==0 && j==0 )
					maxVal[k] = minVal[k] = val;
				if( val > maxVal[k] )
					maxVal[k] = val;
				if( val < minVal[k] )
					minVal[k] = val;
			}

	for(uint32_t k=0;k<n;k++)
		spread[k] = maxVal[k] - minVal[k];

	delete[] maxVal;
	delete[] minVal;

	W.setZero(n,d);

	cout << "Initial Condition: ( ";

	for(uint i=0;i<d;i++) {
		bigVal = 0.0;
		bigAxis = 0;
		for(uint32_t k=0;k<n;k++) {
			if( spread[k] > bigVal ) {
				bigVal = spread[k];
				bigAxis = k;
			}
		}
		spread[bigAxis] = 0.0;
		cout << bigAxis;
		if( i<d-1 )  cout << ", ";
		W(bigAxis,i) = 1.0;
	}
	cout << " )" << endl;

	delete[] spread;

	Wt = W.transpose();
}



double ProjSecantSystem::Find() {
	for(uint32_t i=0;i<data->numFiles;i++)
		if ( !Proj[i]->initialised ) Proj[i]->Init();

	GetInitial();
	
	cost = -1.0;
	minPLen = 0.0;

	//GradientDescent<Grassmannian,MetricFrobenius> OP;
	ConjugateGradient<Grassmannian,MetricFrobenius> OP;
	//QuasiNewton<Grassmannian,MetricFrobenius,LinearMapM> OP;
	MetricFrobenius M;

	OP.maxSteps = maxIts;
	OP.L.obj = this;
	OP.L.metric = &M;
	OP.L.S = ProjSystem::Cost;
	OP.L.gradS = ProjSystem::GradCost;
	OP.L.alpha = 2.0;

	OP.Optimize( W );
	Wt = W.transpose();

	cost = CostBasic(W,minPLen);

	return minPLen;

}


double ProjSecantSystem::CostFunction( const SecantsSystem& secants, const MatrixXd &X ) const {
	double c = 0.0, v = 1.0;
	for(uint32_t i=0;i<data->numFiles;i++) {
		c += ProjSecant::CostFunction(secants);
	}
	c /= data->numFiles;
	return c;
}

MatrixXd ProjSecantSystem::DerivCostBasic( const MatrixXd &Y ) const {
	MatrixXd FY;
	FY.setZero(Y.rows(),Y.cols());
	cout << ",";
	for(uint32_t i=0;i<data->numFiles;i++) {
		FY += Proj[i]->DerivCostBasic(Y);
	}
	FY /= data->numFiles;
	return FY;
}

void ProjSecantSystem::AnalyseSecants() {
	double xMin = 1.0, xMax = 0.0, xMean = 0.0, total = 0.0, len;
	uint32_t n = 0;
	for(uint32_t i=0;i<data->numFiles;i++) {
		for(uint32_t j=0;j<Proj[i]->nSecants;j++) {
			len = ( Wt * Proj[i]->GetSecant(j) ).norm();
			if( len < xMin ) xMin = len;
			if( len > xMax ) xMax = len;
			total += len;
			n++;
		}
	}
	xMean = total / n;
	cout << endl << "Projected Lengths: Range = [ " << xMin << ", " << xMax << " ], Mean = " << xMean << endl;
}

double ProjSecantSystem::Cost( const MatrixXd &X, const void* obj ) {
	const ProjSecantSystem* ps = (const ProjSecantSystem*) obj;
	double m;
	return ps->CostFunction(X,m);
}

MatrixXd ProjSecantSystem::GradCost( const MatrixXd &X, const void* obj ) {
	const ProjSecantSystem* ps = (const ProjSecantSystem*) obj;
	return ProjSecantSystem::ComputeGradient(X,ps->CostFunctionDerivative(X));
}
