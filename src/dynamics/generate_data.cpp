#include <iostream>
#include <fstream>
#include <DRDSP/dynamics/generate_data.h>

using namespace std;
using namespace DRDSP;

DataGenerator::DataGenerator( ModelParameterized& m ) : model(m), rk(model), tStart(0), tInterval(10), print(200), pMin(0), pMax(1), pDelta(0.1), binaryOutput(true), textOutput(false) {
	initial.setZero(m.dimension);
}

void DataGenerator::GenerateFiles() {
	for(double p=pMin;p<=pMax;p+=pDelta) {
		GenerateSingleFile(p);
	}
}

void DataGenerator::GenerateSingleFile( double param ) {
	ofstream outBin, outTxt;
	stringstream fn;
	outTxt.precision(16);
	model.parameter(0) = param;
	fn.str("");
	fn << "data/pp" << param;
	cout << fn.str() << endl;
		
	if( binaryOutput ) {
		fn << ".bin";
		outBin.open(fn.str(),ios::binary);
	}
	if( textOutput ) {
		fn.str("");
		fn << "data/pp" << param << ".csv";
		outTxt.open(fn.str());
	}

	rk.state = initial;
	rk.Advance(tStart);

	double t = 0.0;
	double dtPrint = tInterval / print;
			
	if( binaryOutput ) {
		outBin.write((char*)&print,sizeof(uint32_t));
		outBin.write((char*)&param,sizeof(double));
	}
	if( textOutput )
		outTxt << print << "," << param << endl;

	while( t <= tInterval ) {
			
		if( binaryOutput ) {
			outBin.write((const char*)&rk.state(0),sizeof(double)*rk.state.size());
		}
		if( textOutput ) {
			for(int i=0;i<rk.state.size();i++) {
				outTxt << rk.state(i) << ",";
			}
			outTxt << endl;
		}
		rk.Advance(dtPrint);
		t += dtPrint;
	}
	if( binaryOutput ) outBin.close();
	if( textOutput ) {
		outTxt.close();
	}
}

DataSet DataGenerator::GenerateDataSet( double param ) {

	cout << "Parameter " << param << endl;

	model.parameter(0) = param;
	rk.state = initial;
	rk.Advance(tStart);

	double dtPrint = tInterval / print;

	DataSet data;
	data.Create(print,model.model.dimension);
	for(uint32_t i=0;i<print;i++) {
		data[i] = rk.state;
		rk.Advance(dtPrint);
	}
	return std::move(data);
}

DataSystem DataGenerator::GenerateDataSystem() {
	DataSystem data;
	uint16_t numParms = 1 + (uint16_t)(( pMax - pMin ) / pDelta);
	data.Create(model.model.dimension,numParms,1);
	double p = pMin;
	for(uint16_t i=0;i<numParms;i++) {
		data.parameters[i].setZero(model.model.parameterDimension);
		data.parameters[i](0) = p;
		data.dataSets[i] = GenerateDataSet(p);
		p += pDelta;
	}
	return std::move(data);
}

