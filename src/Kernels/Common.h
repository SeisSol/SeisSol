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

#include "Common/Constants.h"
#include "Kernels/Precision.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include <algorithm>
#include <cassert>
#include <type_traits>
#include <utility>
#include <yateto.h>

/**
 * Uses SFINAE to generate the following functions:
 *
 * has_NAME<T>::value -> true if class T has member NAME and false otherwise
 * set_NAME<T>(kernel, ptr) -> sets kernel.NAME = ptr if class T has member NAME and does nothing
 * otherwise get_static_ptr_NAME<T>() returns &T::NAME[0] if class T has member NAME and nullptr
 * otherwise get_ptr_NAME<T>(T& obj) returns &obj.NAME[0] if class T has member NAME and nullptr
 * otherwise
 */
#define GENERATE_HAS_MEMBER(NAME)                                                                  \
  namespace seissol::kernels {                                                                     \
  template <typename T>                                                                            \
  struct has_##NAME {                                                                              \
    template <typename U>                                                                          \
    static constexpr decltype(std::declval<U>().NAME, bool()) test(int) {                          \
      return true;                                                                                 \
    }                                                                                              \
    template <typename U>                                                                          \
    static constexpr bool test(...) {                                                              \
      return false;                                                                                \
    }                                                                                              \
    static constexpr bool value = test<T>(int());                                                  \
  };                                                                                               \
  template <class T>                                                                               \
  auto set_##NAME(T& kernel, decltype(T::NAME) ptr) ->                                             \
      typename std::enable_if_t<has_##NAME<T>::value> {                                            \
    kernel.NAME = ptr;                                                                             \
  }                                                                                                \
  template <class T, class S>                                                                      \
  auto set_##NAME(T&, const S&) -> typename std::enable_if_t<!has_##NAME<T>::value> {}             \
  template <class T>                                                                               \
  constexpr auto get_static_ptr_##NAME() ->                                                        \
      typename std::enable_if_t<has_##NAME<T>::value, decltype(&T::NAME[0])> {                     \
    return &T::NAME[0];                                                                            \
  }                                                                                                \
  template <class T>                                                                               \
  constexpr auto get_static_ptr_##NAME() ->                                                        \
      typename std::enable_if_t<!has_##NAME<T>::value, void*> {                                    \
    return nullptr;                                                                                \
  }                                                                                                \
  template <class T>                                                                               \
  constexpr auto get_ptr_##NAME(T& obj) ->                                                         \
      typename std::enable_if_t<has_##NAME<T>::value, decltype(&obj.NAME[0])> {                    \
    return &obj.NAME[0];                                                                           \
  }                                                                                                \
  template <class T>                                                                               \
  constexpr auto get_ptr_##NAME(T&) -> typename std::enable_if_t<!has_##NAME<T>::value, void*> {   \
    return nullptr;                                                                                \
  }                                                                                                \
  template <class T>                                                                               \
  constexpr auto get_ref_##NAME(T& obj) ->                                                         \
      typename std::enable_if_t<has_##NAME<T>::value, decltype(obj.NAME)&> {                       \
    return obj.NAME;                                                                               \
  }                                                                                                \
  template <class T>                                                                               \
  constexpr auto get_ref_##NAME(T& obj) -> typename std::enable_if_t<!has_##NAME<T>::value, T&> {  \
    return obj;                                                                                    \
  }                                                                                                \
  template <class T>                                                                               \
  constexpr auto get_cref_##NAME(const T& obj) ->                                                  \
      typename std::enable_if_t<has_##NAME<T>::value, const decltype(obj.NAME)&> {                 \
    return obj.NAME;                                                                               \
  }                                                                                                \
  template <class T>                                                                               \
  constexpr auto get_cref_##NAME(const T& obj) ->                                                  \
      typename std::enable_if_t<!has_##NAME<T>::value, const T&> {                                 \
    return obj;                                                                                    \
  }                                                                                                \
  }

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
  unsigned int numberOfBasisFunctions = getNumberOfBasisFunctions(convergenceOrder);
  return getNumberOfAlignedReals<RealT>(numberOfBasisFunctions);
}

/**
 * Get the # of derivatives of basis functions aligned to the given boundaries.
 *
 * @param convergenceOrder convergence order.
 * @param alignment alignment in bytes.
 * @return aligned number of basis functions.
 **/
template <typename RealT = real>
constexpr unsigned
    getNumberOfAlignedDerivativeBasisFunctions(unsigned int convergenceOrder = ConvergenceOrder,
                                               unsigned int alignment = Vectorsize) {
  return (convergenceOrder > 0)
             ? getNumberOfAlignedBasisFunctions<RealT>(convergenceOrder) +
                   getNumberOfAlignedDerivativeBasisFunctions<RealT>(convergenceOrder - 1)
             : 0;
}

template <typename FamilyT, typename ContainerT>
void setupContainer(ContainerT& container, const real* data) {
  for (std::size_t i = 0; i < yateto::numFamilyMembers<FamilyT>(); ++i) {
    container.data[i] = data + yateto::computeFamilySize<FamilyT>(1, i);
  }
}

/**
 * uses SFINAE to check if class T has a size() function.
 */
template <typename T>
struct HasSize {
  template <typename U>
  static constexpr decltype(U::Size, bool()) test(int /*unused*/) {
    return true;
  }
  template <typename U>
  static constexpr bool test(...) {
    return false;
  }
  static constexpr bool Value = test<T>(int());
};

/**
 * returns T::size() if T has size function and 0 otherwise
 */
template <class T>
constexpr auto size() -> std::size_t {
  if constexpr (HasSize<T>::Value) {
    if constexpr (std::is_array_v<decltype(T::Size)>) {
      return yateto::computeFamilySize<T>();
    } else {
      return T::size();
    }
  } else {
    return 0;
  }
}
} // namespace kernels

constexpr bool isDeviceOn() {
#ifdef ACL_DEVICE
  return true;
#endif
  return false;
}
} // namespace seissol

// for now, make these #defines constexprs. Soon, they should be namespaced.
constexpr std::size_t NumBasisFunctions = seissol::kernels::getNumberOfBasisFunctions();
constexpr std::size_t NumAlignedBasisFunctions =
    seissol::kernels::getNumberOfAlignedBasisFunctions();
constexpr std::size_t NumAlignedDerivativeBasisFunctions =
    seissol::kernels::getNumberOfAlignedDerivativeBasisFunctions();

#endif // SEISSOL_SRC_KERNELS_COMMON_H_
