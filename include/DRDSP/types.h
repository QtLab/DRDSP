#ifndef INCLUDED_TYPES
#define INCLUDED_TYPES

#include <cstdint>

#pragma warning( disable : 4714 ) // function '...' marked as __forceinline not inlined
#pragma warning( disable : 4800 ) // 'int' : forcing value to bool 'true' or 'false' (performance warning)
#pragma warning( disable : 4510 ) // default constructor could not be generated
#pragma warning( disable : 4610 ) // can never be instantiated - user defined constructor required

#include <Eigen/Core>

#pragma warning( default : 4800 )
#pragma warning( default : 4510 )
#pragma warning( default : 4610 )

namespace Eigen {
	typedef Matrix<double,1,1> Vector1d;
}

using namespace Eigen;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#endif
