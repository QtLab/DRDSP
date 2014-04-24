#include <DRDSP/data/histogram.h>
#include <DRDSP/misc.h>
#include <fstream>
#include <numeric>
#include <algorithm>

using namespace DRDSP;

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
	
void Histogram::Create( uint32_t nBins ) {
	bins.resize(nBins);
}
	
void Histogram::Destroy() {
	bins = vector<Bin>();
}
	
uint32_t Histogram::TotalFrequency() const {
	return accumulate(
		begin(bins),
		end(bins),
		uint32_t(0),
		[]( uint32_t x, const Bin& y ){ return x + y.frequency; }
	);
}

void Histogram::WriteCSV( const char* filename ) const {
	ofstream out(filename);
	if( !out ) return;
	uint32_t total = TotalFrequency();
	for(const auto& bin : bins) {
		double freqDensity = bin.FrequencyDensity();
		out << bin.Centre() << ",";
		out << bin.frequency << ",";
		out << (double)bin.frequency / total << ",";
		out << freqDensity << ",";
		out << freqDensity / total << endl;
	}
}

Histogram HistogramGenerator::Generate( const vector<double>& data ) const {
	return Generate( data.data(), data.size() );
}

Histogram HistogramGenerator::Generate( const double* data, size_t N ) const {
	
	uint32_t nBins = numBins;
	if( nBins == 0 ) {
		nBins = uint32_t(1 + N / 20);
	}
		
	Histogram H(nBins);
	
	auto minmax = minmax_element(data,data+N);
	
	double minVal = *minmax.first;
	double maxVal = *minmax.second;

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
			for(auto& bin : H.bins) {
				if( bin.Test( data[i] ) ) {
					++bin;
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
			if( j >= 0 && j < (int32_t)nBins ) {
				++H.bins[j];
			}
		}
	}

	return H;
}


