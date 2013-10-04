#ifndef INCLUDED_TYPES
#define INCLUDED_TYPES

#ifdef _MSC_VER
	typedef __int8 int8_t;
	typedef __int16 int16_t;
	typedef __int32 int32_t;
	typedef __int64 int64_t;
	typedef unsigned __int8 uint8_t;
	typedef unsigned __int16 uint16_t;
	typedef unsigned __int32 uint32_t;
	typedef unsigned __int64 uint64_t;
#else
	#include <stdint.h>
#endif

//#define nullptr 0
typedef unsigned int uint;

#include <Eigen/Core>
using namespace Eigen;

#define M_PI		3.141592653589793238462

#endif
