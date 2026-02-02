// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_CONDITION_H_
#define SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_CONDITION_H_

#include "EncodedConstants.h"

#include <assert.h>
#include <type_traits>

namespace seissol::recording {
template <typename T>
constexpr bool isEncodedConstant() {
  return std::is_same_v<FaceKinds, T> || std::is_same_v<KernelNames, T> ||
         std::is_same_v<FaceId, T> || std::is_same_v<FaceRelations, T> ||
         std::is_same_v<DrFaceRelations, T> || std::is_same_v<ComputationKind, T> ||
         std::is_same_v<ExchangeInfo, T> || std::is_same_v<inner_keys::Wp::Id, T> ||
         std::is_same_v<inner_keys::Dr::Id, T> || std::is_same_v<inner_keys::Indices::Id, T> ||
         std::is_same_v<inner_keys::Material::Id, T>;
}

template <class T, std::enable_if_t<isEncodedConstant<T>()>>
class Condition {
  public:
  Condition() = delete;

  explicit Condition(T initialEncoding)
      : highBitsMask_(~((~static_cast<size_t>(0)) << static_cast<size_t>(T::Count))),
        encoding_(static_cast<size_t>(initialEncoding)) {}

  Condition& operator!() {
    encoding_ = highBitsMask_ & (~encoding_);
    return *this;
  }

  Condition& operator||(const Condition& other) {
    encoding_ = encoding_ | other.encoding_;
    return *this;
  }

  Condition& negate() { return !(*this); }

  size_t getEncoding() { return encoding_; }

  private:
  size_t highBitsMask_;
  size_t encoding_;
};
} // namespace seissol::recording

/** Implements "OR" operation.
 *
 * You won't be able negate the condition after this.
 * Use Condition class to perform more sophisticated logical operations.
 *
 * NOTE: The function was designed to handle simple condition encodings
 * Refer to Condition Class if you need much more sophisticated behaviour
 */
template <typename T>
std::enable_if_t<seissol::recording::isEncodedConstant<T>(), size_t> operator||(const T& lhs,
                                                                                const T& rhs) {
  return (static_cast<size_t>(lhs) | static_cast<size_t>(rhs));
}

/** Returns the actual value of enum item. Behavior is similar as pointer dereference.
 *
 * NOTE: The function was designed to handle simple condition encodings
 * Refer to Condition Class if you need much more sophisticated behaviour
 */
template <typename T>
constexpr std::enable_if_t<seissol::recording::isEncodedConstant<T>(), size_t>
    operator*(const T& condition) {
  return static_cast<size_t>(condition);
}

/** Implements negation operation.
 *
 * NOTE: The function was designed to handle simple condition encodings
 * Refer to Condition Class if you need much more sophisticated behaviour
 */
template <typename T>
std::enable_if_t<seissol::recording::isEncodedConstant<T>(), size_t> operator!(const T& condition) {
  const size_t highBitsMask = ~((~static_cast<size_t>(0)) << static_cast<size_t>(T::Count));
  return highBitsMask & (~static_cast<size_t>(condition));
}

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_CONDITION_H_
