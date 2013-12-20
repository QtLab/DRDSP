#ifndef INCLUDED_DATA_HISTOGRAM
#define INCLUDED_DATA_HISTOGRAM
#include <fstream>
#include "../types.h"

using namespace std;

namespace DRDSP {

	/*!
	 * \brief A single bin for the histogram.
	 */
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

	/*!
	 * \brief A histogram.
	 *
	 * This is just an array of bins with some helper functions.
	 */
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
		uint32_t TotalFrequency() const; //! Sums the bin frequencies
		void WriteCSV( const char* filename ) const;
	};

	/*!
	 * \brief A class for generating histograms from data.
	 */
	struct HistogramGenerator {
		uint32_t numBins; //! The number of bins that we want the histogram to contain
		double clampMin, clampMax; //! clamping limits the range of the bins
		bool clamp,    //! Perform clamping
			 logScale; //! Use a logarithmic scale for the bin sizes
	
		HistogramGenerator();
		HistogramGenerator( uint32_t nBins );
		Histogram Generate( const double* data, size_t N ) const; //! Generates a histogram from the given data
	};

}

#endif
