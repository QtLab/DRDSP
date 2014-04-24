#ifndef INCLUDED_DATA_HISTOGRAM
#define INCLUDED_DATA_HISTOGRAM
#include <stdint.h>
#include <vector>

using namespace std;

namespace DRDSP {

	/*!
	 * \brief A single bin for the histogram.
	 */
	struct Bin {
		double minValue, maxValue;
		uint32_t frequency;

		Bin() : minValue(0.0), maxValue(0.0), frequency(0) {}
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
		vector<Bin> bins;
	
		Histogram() = default;
		Histogram( uint32_t nBins ) : bins(nBins) {}
		void Create( uint32_t nBins );
		void Destroy();
		uint32_t TotalFrequency() const;              //!< Sums the bin frequencies
		void WriteCSV( const char* filename ) const;
	};

	/*!
	 * \brief A class for generating histograms from data.
	 */
	struct HistogramGenerator {
		uint32_t numBins;               //!< The number of bins that we want the histogram to contain
		double clampMin, clampMax;      //!< clamping limits the range of the bins
		bool clamp,                     //!< Perform clamping
			 logScale;                  //!< Use a logarithmic scale for the bin sizes
	
		HistogramGenerator() : numBins(0), clampMin(0.0), clampMax(0.0), clamp(false), logScale(false) {}
		HistogramGenerator( uint32_t nBins ) : numBins(nBins) {}
		Histogram Generate( const double* data, size_t N ) const;   //!< Generates a histogram from the given data
		Histogram Generate( const vector<double>& data ) const;     //!< Generates a histogram from the given data
	};

}

#endif
