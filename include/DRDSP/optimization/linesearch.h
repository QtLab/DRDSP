#ifndef INCLUDED_OPTIMIZATION_LINESEARCH
#define INCLUDED_OPTIMIZATION_LINESEARCH
#include <stdint.h>
#include <iostream>
#include <cmath>
#include "../misc.h"

using namespace std;

namespace DRDSP {

	struct LineSearch {
		double alpha = 1.0e-2,
		       alphaMax = 1.0e6,
		       mu = 0.00001,
		       eta = 0.9;

		double S0, DS0;

		template<typename S,typename DS>
		double Search( const S& cost, const DS& costDeriv ) {
			static const double ratio = 2.0;
			static const uint8_t maxIterations = 5;

			DS0 = costDeriv(0.0);
			
			if( DS0 >= 0.0 ) {
				cout << "LineSearch: Bad Tangent Vector." << endl;
				return 0.0;
			}

			S0 = cost(0.0);
			double a = 0.0;
			double ga = DS0;
			double c = alpha / ratio;
			double gc;
			do {
				c = Clamp( c * ratio, 0.0, alphaMax );
				gc = costDeriv(c);
				if( c == alphaMax ) {
					if( gc < 0.0 ) {
						if( cost(c) < S0 ) {
							alpha = c;
							return alpha;
						} else return 0.0;
					}
					break;
				}
			} while( gc < 0.0 );

			double b = c;
			double gb = gc;
			double Sc = 0.0;
			for(uint8_t i=0;i<maxIterations;++i) {
				c = Secant(a,b,ga,gb);
				gc = costDeriv(c);
				Sc = cost(c);

				if( Wolfe(Sc,gc) ) {
					alpha = c;
					return alpha;
				}
				if( gc < 0.0 ) {
					a = c;
					ga = gc;
				} else if( gc > 0.0 ) {
					b = c;
					gb = gc;
				} else break;
			}
			if( Sc < S0 ) {
				alpha = c;
				return alpha;
			} else {
				return 0.0;
			}
		}

	protected:

		bool Wolfe( double Sc, double gc ) const {

			double c1 = Sc - S0;

			if( c1 > 0.0 ) return false;

			if( DS0 == 0.0 ) return false;

			c1 /= DS0;

			if( c1 > mu ) return false;

			double d2 = std::abs( DS0 );

			if( d2 == 0.0 ) return false;

			double c2 = std::abs( gc ) / d2;

			if( c2 > eta ) return false;

			return true;
		}

		static double Secant( double a, double b, double ga, double gb ) {
			return a - ga*(b-a)/(gb-ga);
		}

	};

}

#endif
