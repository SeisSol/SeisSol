#ifndef SEISSOL_SYCL_FUNCTIONS_H
#define SEISSOL_SYCL_FUNCTIONS_H

#include <array>
#include <cstddef>
#include <cstdint>

namespace seissol::functions {
/**
 * @brief sycl device standard math functions used in template metaprogramming
 */
struct SyclStdFunctions {
  template <typename T>
  static inline T exp(T value) {
    return std::exp(value);
  }
  template <typename T1, typename... T>
  static inline T1 max(T1 value1, T... value) {
    return std::max(value1, value...);
  }
  template <typename T1, typename... T>
  static inline T1 min(T1 value1, T... value) {
    return std::min(value1, value...);
  }
  template <typename T>
  static inline T ceil(T value) {
    return std::ceil(value);
  }
  template <typename T>
  static inline T floor(T value) {
    return std::floor(value);
  }
};
} // namespace seissol::functions

#endif // SEISSOL_SYCL_FUNCTIONS_H
