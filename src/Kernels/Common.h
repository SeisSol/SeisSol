// SPDX-FileCopyrightText: 2014 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_COMMON_H_
#define SEISSOL_SRC_KERNELS_COMMON_H_

#include "Alignment.h"
#include "Common/Constants.h"
#include "Common/Typedefs.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "Kernels/Precision.h"

#include <algorithm>
#include <cassert>
#include <type_traits>
#include <utility>

// TODO: once we use C++20, remove the SFINAE usage.

/**
 * Generate the following functions (needs a macro, since it's qualifier dependent)
 *
 * has_NAME<T>::value -> true if class T has member NAME and false otherwise
 *
 * set_NAME<T>(kernel, ptr) -> sets kernel.NAME = ptr if class T has member NAME and does nothing
 * otherwise
 *
 * get_static_ptr_NAME<T>() returns &T::NAME[0] if class T has member NAME and nullptr otherwise
 *
 * get_ptr_NAME<T>(T& obj) returns &obj.NAME[0] if class T has member NAME and nullptr
 * otherwise
 */
#define GENERATE_HAS_MEMBER(NAME)                                                                  \
  namespace {                                                                                      \
  using namespace seissol::kernels;                                                                \
  template <typename T, typename = void>                                                           \
  struct has_##NAME {                                                                              \
    static constexpr bool value = false;                                                           \
  };                                                                                               \
  template <typename T>                                                                            \
  struct has_##NAME<T, decltype(std::declval<T>().NAME, void())> {                                 \
    static constexpr bool value = true;                                                            \
  };                                                                                               \
  template <typename T, typename PtrT>                                                             \
  void set_##NAME(T& kernel, PtrT&& ptr) {                                                         \
    if constexpr (has_##NAME<T>::value) {                                                          \
      kernel.NAME = std::forward<PtrT>(ptr);                                                       \
    }                                                                                              \
  }                                                                                                \
  template <typename T>                                                                            \
  constexpr auto get_static_ptr_##NAME() {                                                         \
    if constexpr (has_##NAME<T>::value) {                                                          \
      return &T::NAME[0];                                                                          \
    } else {                                                                                       \
      return nullptr;                                                                              \
    }                                                                                              \
  }                                                                                                \
  template <typename T>                                                                            \
  constexpr auto get_ptr_##NAME(T& obj) {                                                          \
    if constexpr (has_##NAME<T>::value) {                                                          \
      return &obj.NAME[0];                                                                         \
    } else {                                                                                       \
      return nullptr;                                                                              \
    }                                                                                              \
  }                                                                                                \
  } // namespace

namespace seissol {
namespace kernels {
/**
 * Gets the number of basis functions for the given convergence order.
 *
 * @param convergenceOrder convergence order.
 * @return number of basis funcitons.
 **/
constexpr unsigned int getNumberOfBasisFunctions(unsigned int convergenceOrder = ConvergenceOrder) {
  return convergenceOrder * (convergenceOrder + 1) * (convergenceOrder + 2) / 6;
}

/**
 * Gets the number of aligned reals, i.e. the number padded to the size of the alignment.
 *
 * @param alignment alignment in bytes.
 * @return aligned number of reals.
 **/
template <typename RealT = real>
constexpr unsigned int getNumberOfAlignedReals(unsigned int numberOfReals,
                                               unsigned int alignment = Vectorsize) {
  // in principle, we could simplify this formula by substituting alignment = alignment /
  // sizeof(real). However, this will cause errors, if alignment is not dividable by sizeof(real)
  // which could happen e.g. if alignment < sizeof(real), or if we have real == long double (if
  // there is ever such a use case, and if the alignment then still makes much sense).
  return (numberOfReals * sizeof(RealT) +
          (alignment - (numberOfReals * sizeof(RealT)) % alignment) % alignment) /
         sizeof(RealT);
}

/**
 * Get the # of basis functions aligned to the given boundaries.
 *
 * @param convergenceOrder convergence order.
 * @param alignment alignment in bytes.
 * @return aligned number of basis functions.
 **/
template <typename RealT = real>
constexpr unsigned int
    getNumberOfAlignedBasisFunctions(unsigned int convergenceOrder = ConvergenceOrder,
                                     unsigned int alignment = Vectorsize) {
  // return (numberOfBasisFunctions(O) * REAL_BYTES + (ALIGNMENT - (numberOfBasisFunctions(O) *
  // REAL_BYTES) % ALIGNMENT) % ALIGNMENT) / REAL_BYTES
  const auto numberOfBasisFunctions = getNumberOfBasisFunctions(convergenceOrder);
  return getNumberOfAlignedReals<RealT>(numberOfBasisFunctions);
}

/**
 * Get the # of derivatives of basis functions aligned to the given boundaries.
 *
 * @param convergenceOrder convergence order.
 * @param alignment alignment in bytes.
 * @return aligned number of basis functions.
 **/
constexpr unsigned
    getNumberOfAlignedDerivativeBasisFunctions(unsigned int convergenceOrder = ConvergenceOrder,
                                               unsigned int alignment = Vectorsize) {
  return (convergenceOrder > 0)
             ? getNumberOfAlignedBasisFunctions(convergenceOrder) +
                   getNumberOfAlignedDerivativeBasisFunctions(convergenceOrder - 1)
             : 0;
}

/**
 * Check if a type has a .size() member.
 */
template <typename T, typename = void>
struct HasSize {
  static constexpr bool Value = false;
  using Type = std::size_t;
};

template <typename T>
struct HasSize<T, decltype(std::declval<T>().size(), void())> {
  static constexpr bool Value = true;
  using Type = decltype(std::declval<T>().size());
};

/**
 * returns T::size() if T has size function and 0 otherwise.
 */
template <class T>
constexpr auto size() -> typename HasSize<T>::Type {
  if constexpr (HasSize<T>::Value) {
    return T::size();
  } else {
    return static_cast<typename HasSize<T>::Type>(0);
  }
}

} // namespace kernels

constexpr bool isDeviceOn() { return HardwareSupport == BuildType::Gpu; }
} // namespace seissol

// for now, make these #defines constexprs. Soon, they should be namespaced.
constexpr std::size_t NumBasisFunctions = seissol::kernels::getNumberOfBasisFunctions();
constexpr std::size_t NumAlignedBasisFunctions =
    seissol::kernels::getNumberOfAlignedBasisFunctions();
constexpr std::size_t NumAlignedDerivativeBasisFunctions =
    seissol::kernels::getNumberOfAlignedDerivativeBasisFunctions();

#endif // SEISSOL_SRC_KERNELS_COMMON_H_
