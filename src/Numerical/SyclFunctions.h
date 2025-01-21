// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_SYCLFUNCTIONS_H_
#define SEISSOL_SRC_NUMERICAL_SYCLFUNCTIONS_H_

#include <CL/sycl.hpp>
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
    return cl::sycl::exp(value);
  }
  template <typename T1, typename... T>
  static inline T1 max(T1 value1, T... value) {
    return cl::sycl::max(value1, value...);
  }
  template <typename T1, typename... T>
  static inline T1 min(T1 value1, T... value) {
    return cl::sycl::min(value1, value...);
  }
  template <typename T>
  static inline T ceil(T value) {
    return cl::sycl::ceil(value);
  }
  template <typename T>
  static inline T floor(T value) {
    return cl::sycl::floor(value);
  }
};
} // namespace seissol::functions

#endif // SEISSOL_SRC_NUMERICAL_SYCLFUNCTIONS_H_
