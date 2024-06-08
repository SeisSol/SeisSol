
#ifndef PRE_PROCESSOR_MACROS_HPP
#define PRE_PROCESSOR_MACROS_HPP

#include "Monitoring/instrumentation.hpp"
#include <cstddef>

#include "Common/constants.hpp"

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

// for now, turn all #defines which are left into constexprs

// workaround for old NVHPC versions (the output would cause errors there)
#ifdef __NVCOMPILER
// we'll leave the comment in the next line in for now, until a NVHPC version is fixed
#if 0 // (__NVCOMPILER_MAJOR__ > 24) || (__NVCOMPILER_MAJOR__ == 24 && __NVCOMPILER_MINOR__ >= 3)
#define NVHPC_AVOID_OMP 0
#else
#define NVHPC_AVOID_OMP 1
#endif
#else
#define NVHPC_AVOID_OMP 0
#endif

#endif
