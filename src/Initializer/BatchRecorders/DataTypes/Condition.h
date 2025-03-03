// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
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

namespace seissol::initializer::recording {
template <typename T>
constexpr bool isEncodedConstant() {
  return std::is_same_v<FaceKinds, T> || std::is_same_v<KernelNames, T> ||
         std::is_same_v<FaceId, T> || std::is_same_v<FaceRelations, T> ||
         std::is_same_v<DrFaceRelations, T> || std::is_same_v<ComputationKind, T> ||
         std::is_same_v<ExchangeInfo, T> || std::is_same_v<inner_keys::Wp::Id, T> ||
         std::is_same_v<inner_keys::Dr::Id, T> || std::is_same_v<inner_keys::Indices::Id, T> ||
         std::is_same_v<inner_keys::Material::Id, T>;
}

template <class T, typename std::enable_if<isEncodedConstant<T>()>::type>
class Condition {
  public:
  Condition() = delete;

  Condition(T initialEncoding) : encoding(static_cast<size_t>(initialEncoding)) {
    highBitsMask = ~((~size_t(0)) << static_cast<size_t>(T::Count));
  }

  Condition& operator!() {
    encoding = highBitsMask & (~encoding);
    return *this;
  }

  Condition& operator||(const Condition& other) {
    encoding = encoding | other.encoding;
    return *this;
  }

  Condition& negate() { return !(*this); }

  size_t getEncoding() { return encoding; }

  private:
  size_t highBitsMask;
  size_t encoding;
  size_t count;
};
} // namespace seissol::initializer::recording

using namespace seissol::initializer::recording;
using namespace seissol;

/** Implements "OR" operation.
 *
 * You won't be able negate the condition after this.
 * Use Condition class to perform more sophisticated logical operations.
 *
 * NOTE: The function was designed to handle simple condition encodings
 * Refer to Condition Class if you need much more sophisticated behaviour
 */
template <typename T>
typename std::enable_if<isEncodedConstant<T>(), size_t>::type operator||(const T& lhs,
                                                                         const T& rhs) {
  return (static_cast<size_t>(lhs) | static_cast<size_t>(rhs));
}

/** Returns the actual value of enum item. Behavior is similar as pointer dereference.
 *
 * NOTE: The function was designed to handle simple condition encodings
 * Refer to Condition Class if you need much more sophisticated behaviour
 */
template <typename T>
constexpr typename std::enable_if<isEncodedConstant<T>(), size_t>::type
    operator*(const T& condition) {
  return static_cast<size_t>(condition);
}

/** Implements negation operation.
 *
 * NOTE: The function was designed to handle simple condition encodings
 * Refer to Condition Class if you need much more sophisticated behaviour
 */
template <typename T>
typename std::enable_if<isEncodedConstant<T>(), size_t>::type operator!(const T& condition) {
  size_t highBitsMask = ~((~size_t(0)) << static_cast<size_t>(T::Count));
  return highBitsMask & (~static_cast<size_t>(condition));
}

#endif // SEISSOL_SRC_INITIALIZER_BATCHRECORDERS_DATATYPES_CONDITION_H_
