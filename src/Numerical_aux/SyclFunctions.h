#ifndef SEISSOL_SYCL_FUNCTIONS_H
#define SEISSOL_SYCL_FUNCTIONS_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <CL/sycl.hpp>

namespace seissol::functions {
/**
 * @brief sycl device standard math functions used in template metaprogramming
 */
struct SyclStdFunctions {
  template <typename T>
  static T exp(T value) { return cl::sycl::exp(value); }
};
}

#endif // SEISSOL_SYCL_FUNCTIONS_H