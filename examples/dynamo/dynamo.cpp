#include "dynamo.h"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace DRDSP;

// Precompute data that doesn't depend on the state or parameter.
void Dynamo::Init() {
	sinheta.setZero(nI);
	cosheta.setZero(nI);
	cotheta.setZero(nI);
	sintheta.setZero(nJ);
	costheta.setZero(nJ);
	c.setZero(nI,nJ);
	pi32.setZero(nI,nJ);
	Bound1.setZero(nJ,nJ);
	Bound2.setZero(nJ,nJ);
	trig1.setZero(nI,nJ);

	jp1.setZero(nJ);
	jm1.setZero(nJ);

	double ebd = ::exp(eta0), temp;
	cosheta(0) = 0.0;
	sinheta(0) = 0.0;
	cotheta(0) = 1.0;
	for(uint32_t i=1;i<nI;i++) {
        temp = s(i);
        cosheta(i) = 0.5*(ebd/temp + temp/ebd);
        sinheta(i) = 0.5*(ebd/temp - temp/ebd);
		cotheta(i) = cosheta(i)/sinheta(i);
	}

	costheta(0) = 1.0;
	sintheta(0) = 0.0;
	for(uint32_t j=1;j<nJ;j++) {
		temp = theta(j);
		costheta(j) = ::cos(temp);
		sintheta(j) = ::sin(temp);
	}

	for(uint32_t j=0;j<nJ;j++) {
		jp1(j) = int(j)+1;
		jm1(j) = int(j)-1;
	}
	jp1(nJ-1) = 0;
	jm1(0) = nJ-1;
	
	for(uint32_t j=0;j<nJ;j++)
		for(uint32_t i=0;i<nI;i++)
			c(i,j) = cosheta(i) - costheta(j);

	
	for(uint32_t j=0;j<nJ;j++) {
		trig1(0,j) = -costheta(j);
		pi32(0,j) = 1.0;
		for(uint32_t i=1;i<nI;i++) {
			temp = c(i,j) / sinheta(i);
			pi32(i,j) = ::sqrt( temp*temp*temp );
			trig1(i,j) = (1.0-cosheta(i)*costheta(j)) / sinheta(i);
		}
	}
	
	ifstream in("jepps.txt");
	if( !in ) cout << "ERROR: jepps.txt fail." << endl;
	for(uint32_t k=0;k<nJ;k++)
		for(uint32_t j=0;j<nJ;j++)
			in >> Bound2(j,k);

	for(uint32_t k=0;k<nJ;k++)
		for(uint32_t j=0;j<nJ;j++)
			in >> Bound1(j,k);

}
