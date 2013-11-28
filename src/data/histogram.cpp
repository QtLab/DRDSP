#include <DRDSP/data/histogram.h>
#include <DRDSP/misc.h>
#include <cmath>

using namespace DRDSP;

Bin::Bin() : minValue(0.0), maxValue(0.0), frequency(0) {}

double Bin::Width() const {
	return maxValue - minValue;
}
	
double Bin::Centre() const {
	return 0.5*(minValue+maxValue);
}
	
double Bin::FrequencyDensity() const {
	return double(frequency) / Width();
}
	
bool Bin::Test( double x ) const {
	return ( x >= minValue ) && ( x < maxValue );
}
	
Bin& Bin::operator++() {
	++frequency;
	return *this;
}
	
Bin& Bin::operator--() {
	--frequency;
	return *this;
}
	
Histogram::Histogram() : bins(nullptr), numBins(0) {}
	
Histogram::Histogram( uint32_t nBins ) : bins(nullptr), numBins(0) {
	Create(nBins);
}
	
Histogram::Histogram( const Histogram& rhs ) {
	Create(rhs.numBins);
	for(uint32_t i=0;i<numBins;i++)
		bins[i] = rhs.bins[i];
}
	
Histogram::Histogram( Histogram&& rhs ) : numBins(rhs.numBins), bins(rhs.bins) {
	rhs.numBins = 0;
	rhs.bins = nullptr;
}
	
Histogram::~Histogram() {
	Destroy();
}
	
Histogram& Histogram::operator=( const Histogram& rhs ) {
	Destroy();
	Create(rhs.numBins);
	for(uint32_t i=0;i<numBins;i++)
		bins[i] = rhs.bins[i];
	return *this;
}
	
Histogram& Histogram::operator=( Histogram&& rhs ) {
	if( this != &rhs ) {
		Destroy();
		numBins = rhs.numBins;
		bins = rhs.bins;
		rhs.numBins = 0;
		rhs.bins = nullptr;
	}
	return *this;
}
	
void Histogram::Create( uint32_t nBins ) {
	bins = new Bin [nBins];
	numBins = nBins;
}
	
void Histogram::Destroy() {
	delete[] bins;
	bins = nullptr;
	numBins = 0;
}
	
uint32_t Histogram::TotalFrequency() const {
	uint32_t total = 0;
	for(uint32_t i=0;i<numBins;i++)
		total += bins[i].frequency;
	return total;
}
	
void Histogram::WriteCSV( const char* filename ) const {
	ofstream out(filename);
	if( !out ) return;
	uint32_t total = TotalFrequency();
	for(uint32_t i=0;i<numBins;i++) {
		double freqDensity = bins[i].FrequencyDensity();
		out << bins[i].Centre() << ",";
		out << bins[i].frequency << ",";
		out << (double)bins[i].frequency / total << ",";
		out << freqDensity << ",";
		out << freqDensity / total << endl;
	}
}

HistogramGenerator::HistogramGenerator() : numBins(0), clampMin(0.0), clampMax(0.0), clamp(false), logScale(false) {}

HistogramGenerator::HistogramGenerator( uint32_t nBins ) : numBins(nBins) {}

Histogram HistogramGenerator::Generate( const double* data, size_t N ) const {
	
	uint32_t nBins = numBins;
	if( nBins == 0 ) {
		nBins = 1 + N / 20;
	}
		
	Histogram H(nBins);
	double minVal = data[0];
	double maxVal = data[0];
	for(size_t i=1;i<N;i++) {
		if( data[i] > maxVal ) {
			maxVal = data[i];
		}
		if( data[i] < minVal ) {
			minVal = data[i];
		}
	}
	
	if( clamp ) {
		minVal = Clamp(minVal,clampMin,clampMax);
		maxVal = Clamp(maxVal,clampMin,clampMax);
	}

	if( logScale ) {
		double factor = pow( maxVal/minVal, 1.0/nBins );
		H.bins[0].minValue = minVal;
		for(uint32_t i=0;i<nBins-1;i++) {
			H.bins[i].maxValue = H.bins[i].minValue * factor;
			H.bins[i+1].minValue = H.bins[i].maxValue;
		}
		H.bins[nBins-1].maxValue = H.bins[nBins-1].minValue * factor;
		for(size_t i=0;i<N;i++) {
			for(uint32_t j=0;j<nBins;j++) {
				if( H.bins[j].Test( data[i] ) ) {
					++H.bins[j];
					break;
				}
			}
		}
	} else {
		double width = ( maxVal - minVal ) / nBins;
		for(uint32_t i=0;i<nBins;i++) {
			H.bins[i].minValue = minVal + width * i;
			H.bins[i].maxValue = minVal + width * (i+1);
		}
		for(size_t i=0;i<N;i++) {
			double temp = ( data[i] - minVal ) / width;
			int32_t j = (int32_t)temp;
			if( j >= 0 && j < nBins ) {
				++H.bins[j];
			}
		}
	}


	return std::move(H);
}


