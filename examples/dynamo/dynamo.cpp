#pragma warning ( disable : 4503 )
#include "dynamo.h"
#include <iostream>
#include <fstream>

using namespace std;
using namespace DRDSP;

void Dynamo::Init() {
	sinheta.setZero(nI);
	cosheta.setZero(nI);
	cotheta.setZero(nI);
	s.setZero(nI);
	sintheta.setZero(nJ);
	costheta.setZero(nJ);
	c.setZero(nI,nJ);
	pi32.setZero(nI,nJ);
	Bound1.setZero(nJ,nJ);
	Bound2.setZero(nJ,nJ);
	trig1.setZero(nI,nJ);

	jp1.setZero(nJ);
	jm1.setZero(nJ);

	double ebd = ::exp(eta0);
	cosheta(0) = numeric_limits<double>::infinity();
	sinheta(0) = numeric_limits<double>::infinity();
	cotheta(0) = 1.0;
	for(uint32_t i=1;i<nI;++i) {
		double temp = i * ds;
		s(i) = temp;
        cosheta(i) = 0.5*(ebd/temp + temp/ebd);
        sinheta(i) = 0.5*(ebd/temp - temp/ebd);
		cotheta(i) = cosheta(i)/sinheta(i);
		//cout << cosheta(i) << ", " << sinheta(i) << ", " << cotheta(i) << endl;
	}

	costheta(0) = 1.0;
	sintheta(0) = 0.0;
	for(uint32_t j=1;j<nJ;++j) {
		double temp = j * dth;
		costheta(j) = ::cos(temp);
		sintheta(j) = ::sin(temp);
	}

	for(uint32_t j=0;j<nJ;++j) {
		jp1(j) = int(j)+1;
		jm1(j) = int(j)-1;
	}
	jp1(nJ-1) = 0;
	jm1(0) = nJ-1;
	
	for(uint32_t j=0;j<nJ;++j) {
		c(0,j) = cosheta(0);
		for(uint32_t i=1;i<nI;++i)
			c(i,j) = cosheta(i) - costheta(j);
	}

	
	for(uint32_t j=0;j<nJ;++j) {
		trig1(0,j) = -costheta(j);
		pi32(0,j) = 1.0;
		for(uint32_t i=1;i<nI;++i) {
			double temp = c(i,j) / sinheta(i);
			pi32(i,j) = ::sqrt( temp*temp*temp );
			trig1(i,j) = (1.0-cosheta(i)*costheta(j)) / sinheta(i);
		}
	}
	
	ifstream in("jepps.txt");
	if( !in ) { 
		cerr << "Dynamo::Init : jepps.txt fail." << endl;
		return;
	}
	for(uint32_t k=0;k<nJ;++k)
		for(uint32_t j=0;j<nJ;++j)
			in >> Bound2(j,k);

	for(uint32_t k=0;k<nJ;++k)
		for(uint32_t j=0;j<nJ;++j)
			in >> Bound1(j,k);

}

VectorXd Dynamo::InitialCondition() const {
	VectorXd state(dimension);

	for(uint32_t i=0;i<nI;++i) {
		double x = s(i);
		for(uint32_t j=0;j<nJ;++j) {
			state[i*nJ+j] = x*(1.0-x)*(sintheta(j) + costheta(j)) * 1.0e-3;
		}
	}
	double eps = 0.0;
	for(uint32_t i=0;i<nI;++i) {
		double x = s(i);
		for(uint32_t j=0;j<nJ;++j) {
			state[N+i*nJ+j] = x*(1.0-x)*(sintheta(j) + costheta(j)) * 1.0e-3 + eps * x * x * (1.0 - x) * sintheta(j) * costheta(j);
		}
	}
	return state;
}
