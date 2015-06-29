#include <initialization/precision.h>

#ifdef SINGLE_PRECISION
#define REAL_BYTES sizeof(float)
#endif
#ifdef DOUBLE_PRECISION
#define REAL_BYTES sizeof(double)
#endif

#include <Kernels/equations.hpp>
