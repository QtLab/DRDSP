#ifndef INCLUDED_TYPES
#define INCLUDED_TYPES

#include <stdint.h>

typedef unsigned int uint;

#pragma warning( disable : 4714 ) // function '...' marked as __forceinline not inlined

#include <Eigen/Core>
using namespace Eigen;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#endif
