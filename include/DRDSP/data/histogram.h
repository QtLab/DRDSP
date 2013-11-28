#ifndef INCLUDED_DATA_HISTOGRAM
#define INCLUDED_DATA_HISTOGRAM
#include <fstream>
#include "../types.h"

using namespace std;

namespace DRDSP {
	struct Bin {
		double minValue, maxValue;
		uint32_t frequency;

		Bin();
		double Width() const;
		double Centre() const;
		double FrequencyDensity() const;
		bool Test( double x ) const;
		Bin& operator++();
		Bin& operator--();
	};

	struct Histogram {
		Bin* bins;
		uint32_t numBins;
	
		Histogram();
		Histogram( uint32_t nBins );
		Histogram( const Histogram& rhs );
		Histogram( Histogram&& rhs );
		~Histogram();
		Histogram& operator=( const Histogram& rhs );
		Histogram& operator=( Histogram&& rhs );
		void Create( uint32_t nBins );
		void Destroy();
		uint32_t TotalFrequency() const;
		void WriteCSV( const char* filename ) const;
	};

	struct HistogramGenerator {
		uint32_t numBins;
		double clampMin, clampMax;
		bool clamp, logScale;
	
		HistogramGenerator();
		HistogramGenerator( uint32_t nBins );
		Histogram Generate( const double* data, size_t N ) const;
	};

}

#endif
