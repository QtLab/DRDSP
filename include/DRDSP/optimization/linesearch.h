#ifndef INCLUDED_LINESEARCH
#define INCLUDED_LINESEARCH
#include "../types.h"
#include "../Geometry/geodesic.h"
#include "../Geometry/metric.h"

namespace DRDSP {
	template<typename TMet>
	struct LineSearch {
		typedef typename TMet::TVec TVec;
		double (*S)( const TVec &X, const void *obj );
		TVec (*gradS)( const TVec &X, const void *obj );
		const void *obj;
		const TMet *metric;
		TVec gradSx;
		double Sx, alpha;

		LineSearch( Geodesic<TVec> *geodesic ) : G(geodesic), alpha(1e-2), alphaMax(1e6), mu(0.00001), eta(0.9), metric(nullptr), obj(nullptr) {}

		double DerivS( const TVec &X, const TVec &V ) {
			return (*metric)(X)(gradS(X,obj),V);
		}

		double Search() {
			double a,b,c,ga,gb,gc,Sc;
			double ratio = 2.0;

			Sx = S(G->GetPoint(),obj);
			derivSx = (*metric)(G->GetPoint())(gradSx,G->GetTangent()); 
			
			if( derivSx < 0.0 ) {
				a = 0.0;
				ga = derivSx;
			} else {
				cout << "LineSearch: Bad Tangent Vector." << endl;
				return 0.0;
			}
			c = alpha/ratio;

			TVec Xc;

			do {
				c *= ratio;
				if( c >= alphaMax ) {
					c = alphaMax;
					cout << "LineSearch: alphaMax hit." << endl;
				}
				Xc = G->Evaluate(c);
				gc = DerivS(Xc,G->ParallelTranslate(G->GetTangent(),c));
				if( c == alphaMax ) {
					if( gc < 0.0 ) {
						if( S(Xc,obj) < Sx ) {
							alpha = c;
							return alpha;
						} else return 0.0;
					}
					break;
				}
				//cout << "! ";
			} while( gc < 0.0 );

			b = c;
			gb = gc;
			for(uint32_t i=0;i<5;i++) {
				c = Secant(a,b,ga,gb);
				Xc = G->Evaluate(c);
				gc = DerivS(Xc,G->ParallelTranslate(G->GetTangent(),c));
				Sc = S(Xc,obj);
				//cout << "\t " << gc << endl;
				if( Wolfe(Sc,gc) ) {
					//cout << "\t Wolftastic!" << endl;
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
			if( Sc < Sx ) {
				alpha = c;
				return alpha;
			} else {
				//cout << "LineSearch: Fail." << endl;
				return 0.0;
			}
		}

	protected:
		Geodesic<TVec> *G;
		double derivSx, mu, eta, alphaMax;

		bool Wolfe( double Sc, double gc ) {
			double c1, c2, d1, d2;

			c1 = Sc - Sx;
			d1 = derivSx;

			if( c1 > 0.0 )
				return false; // herp derp

			if( d1 != 0.0 )
				c1 /= d1;
			else
				return false;
			//cout << "c1 = " << c1 << endl;
			if( c1 > mu )
				return false;

			c2 = abs( gc );
			d2 = abs( derivSx );

			if( d2 != 0.0 )
				c2 /= d2;
			else
				return false;
			
			//cout << "c2 = " << c2 << endl;
			if( c2 > eta )
				return false;

			return true;
		}

		double Secant( double a, double b, double ga, double gb ) {
			return a - ga*(b-a)/(gb-ga);
		}

		double Riddler( double x1, double x2, double x3, double gx1, double gx2, double gx3 ) {
			double m = ( gx1 >= gx2 )?1.0:-1.0;
			double s = gx3*gx3 - gx1*gx2;
			if( s <= 0.0 ) return x3;
			return x3 + (x3 - x1) * ( m * gx3 )/sqrt( s );
		}

	};
}

#endif

