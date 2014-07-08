#ifndef INCLUDED_TYPES
#define INCLUDED_TYPES

#include <cstdint>

#pragma warning( disable : 4714 ) // function '...' marked as __forceinline not inlined

#include <Eigen/Core>

namespace Eigen {
	typedef Matrix<double,1,1> Vector1d;
}

using namespace Eigen;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#endif
