
#ifndef PRE_PROCESSOR_MACROS_HPP
#define PRE_PROCESSOR_MACROS_HPP

#include "Monitoring/instrumentation.hpp"
#include <cstddef>

#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)

// for now, turn all #defines which are left into constexprs

constexpr std::size_t PAGESIZE_HEAP = 2097152;
constexpr std::size_t PAGESIZE_STACK = 4096;

constexpr std::size_t ALLOW_POSSILBE_ZERO_LENGTH_ARRAY(std::size_t X) {
    return X == 0 ? 1 : X;
}

#endif
