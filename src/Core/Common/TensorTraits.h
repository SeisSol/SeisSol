// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_COMMON_TENSORTRAITS_H_
#define SEISSOL_SRC_CORE_COMMON_TENSORTRAITS_H_

#include <cstddef>
#include <utility>
#include <yateto.h>

/**
 * The traits of the generated tensors: their sizes, and those of their families. A tensor that a
 * configuration lacks is void in its generated code (as the code generator of yateto for several
 * configurations writes it); it has none of the traits, so its size is zero. This is the one place
 * that knows that convention.
 */
namespace seissol::kernels {

/**
 * Check if a type has a .size(...) member; Type is what it returns, or std::size_t without one.
 */
template <typename T, typename... Args>
struct HasSize {
  static constexpr bool Value = false;
  using Type = std::size_t;
};

template <typename T, typename... Args>
  requires requires { std::declval<T>().size(std::declval<Args>()...); }
struct HasSize<T, Args...> {
  static constexpr bool Value = true;
  using Type = decltype(std::declval<T>().size(std::declval<Args>()...));
};

/**
 * returns T::size() if T has size function and 0 otherwise.
 */
template <class T, typename... Args>
constexpr auto size(Args... args) -> HasSize<T, Args...>::Type {
  if constexpr (HasSize<T, Args...>::Value) {
    return T::size(args...);
  } else {
    return static_cast<HasSize<T, Args...>::Type>(0);
  }
}

// A generated family is recognised by its Size array, which is what yateto
// reads to add the members up. Its Container is a class template, so naming it
// without arguments never matched and every family measured as empty.
template <typename T>
struct IsFamilyInternal {
  static constexpr bool Value = requires { T::Size; };
};

template <class T>
constexpr auto familyMembers() -> std::size_t {
  if constexpr (IsFamilyInternal<T>::Value) {
    return yateto::numFamilyMembers<T>();
  } else {
    return static_cast<std::size_t>(0);
  }
}

template <class T>
constexpr auto familySize(std::size_t alignedReals = 1) -> std::size_t {
  if constexpr (IsFamilyInternal<T>::Value) {
    return yateto::computeFamilySize<T>(alignedReals);
  } else {
    return static_cast<std::size_t>(0);
  }
}

} // namespace seissol::kernels

#endif // SEISSOL_SRC_CORE_COMMON_TENSORTRAITS_H_
