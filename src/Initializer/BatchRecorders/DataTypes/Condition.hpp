#ifndef SEISSOL_CONDITION_HPP
#define SEISSOL_CONDITION_HPP

#include "EncodedConstants.hpp"
#include <assert.h>
#include <type_traits>

namespace seissol {
namespace initializers {
namespace recording {

template <typename T> constexpr bool isEncodedConstant() {
  return std::is_same<FaceKinds, T>::value || std::is_same<KernelNames, T>::value ||
         std::is_same<FaceId, T>::value || std::is_same<FaceRelations, T>::value ||
         std::is_same<DrFaceRelations, T>::value || std::is_same<ComputationKind, T>::value ||
         std::is_same<EntityId, T>::value || std::is_same<ExchangeInfo, T>::value;
}


template <class T, typename std::enable_if<isEncodedConstant<T>()>::type> class Condition {
public:
  Condition() = delete;

  Condition(T initialEncoding) : encoding(static_cast<size_t>(initialEncoding)) {
    highBitsMask = ~((~size_t(0)) << static_cast<size_t>(T::Count));
  }

  Condition &operator!() {
    encoding = highBitsMask & (~encoding);
    return *this;
  }

  Condition &operator||(const Condition &other) {
    encoding = encoding | other.encoding;
    return *this;
  }

  Condition &negate() {
    return !(*this);
  }

  size_t getEncoding() {
    return encoding;
  }

private:
  size_t highBitsMask;
  size_t encoding;
  size_t count;
};
} // namespace recording
} // namespace initializers
} // namespace seissol


using namespace seissol::initializers::recording;
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
typename std::enable_if<isEncodedConstant<T>(), size_t>::type operator||(const T &lhs,
                                                                         const T &rhs) {
  return (static_cast<size_t>(lhs) | static_cast<size_t>(rhs));
}


/** Returns the actual value of enum item. Behavior is similar as pointer dereference.
 *
 * NOTE: The function was designed to handle simple condition encodings
 * Refer to Condition Class if you need much more sophisticated behaviour
 */
template <typename T>
constexpr typename std::enable_if<isEncodedConstant<T>(), size_t>::type
operator*(const T &condition) {
  return static_cast<size_t>(condition);
}


/** Implements negation operation.
 *
 * NOTE: The function was designed to handle simple condition encodings
 * Refer to Condition Class if you need much more sophisticated behaviour
 */
template <typename T>
typename std::enable_if<isEncodedConstant<T>(), size_t>::type operator!(const T &condition) {
  size_t highBitsMask = ~((~size_t (0)) << static_cast<size_t>(T::Count));
  return highBitsMask & (~static_cast<size_t>(condition));
}

#endif // SEISSOL_CONDITION_HPP
