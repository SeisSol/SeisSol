// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_TESTHELPER_H_
#define SEISSOL_TESTS_TESTHELPER_H_

#include <doctest.h>

#include "Common/Filesystem.h"
#include "Setup.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <ostream>
#include <type_traits>
#include <vector>

namespace seissol::unit_test {

/// Whether two floating-point values have the same bits: for a result that is to be reproduced
/// exactly, down to the sign of a zero.
template <typename T>
inline bool bitwiseEqual(T lhs, T rhs) {
  static_assert(std::is_floating_point_v<T> && (sizeof(T) == 4 || sizeof(T) == 8),
                "bitwiseEqual compares float and double");
  using Bits = std::conditional_t<sizeof(T) == 8, std::uint64_t, std::uint32_t>;
  Bits lhsBits = 0;
  Bits rhsBits = 0;
  std::memcpy(&lhsBits, &lhs, sizeof(T));
  std::memcpy(&rhsBits, &rhs, sizeof(T));
  return lhsBits == rhsBits;
}

/// The same for `count` values each.
template <typename T>
inline bool bitwiseEqual(const T* lhs, const T* rhs, std::size_t count) {
  for (std::size_t i = 0; i < count; ++i) {
    if (!bitwiseEqual(lhs[i], rhs[i])) {
      return false;
    }
  }
  return true;
}

// Inspired by doctest's Approx, slightly modified
class AbsApprox {
  public:
  inline explicit AbsApprox(double value);

  inline AbsApprox operator()(double value) const;

  inline AbsApprox& epsilon(double newEpsilon);

  inline AbsApprox& delta(double newDelta);

  friend bool operator==(double lhs, const AbsApprox& rhs);

  friend bool operator==(const AbsApprox& lhs, double rhs);

  friend bool operator!=(double lhs, const AbsApprox& rhs);

  friend bool operator!=(const AbsApprox& lhs, double rhs);

  friend bool operator<=(double lhs, const AbsApprox& rhs);

  friend bool operator<=(const AbsApprox& lhs, double rhs);

  friend bool operator>=(double lhs, const AbsApprox& rhs);

  friend bool operator>=(const AbsApprox& lhs, double rhs);

  friend bool operator<(double lhs, const AbsApprox& rhs);

  friend bool operator<(const AbsApprox& lhs, double rhs);

  friend bool operator>(double lhs, const AbsApprox& rhs);

  friend bool operator>(const AbsApprox& lhs, double rhs);

  friend doctest::String toString(const AbsApprox& in);

  private:
  double epsilon_{std::numeric_limits<double>::epsilon()};
  double delta_{0.0};
  double value_{0.0};
};

AbsApprox::AbsApprox(double value) : value_(value) {}

AbsApprox AbsApprox::operator()(double newValue) const {
  AbsApprox approx(newValue);
  approx.epsilon(epsilon_);
  approx.delta(delta_);
  return approx;
}

AbsApprox& AbsApprox::epsilon(double newEpsilon) {
  epsilon_ = newEpsilon;
  return *this;
}

AbsApprox& AbsApprox::delta(double newDelta) {
  delta_ = newDelta;
  return *this;
}

inline bool operator==(double lhs, const AbsApprox& rhs) {
  return std::abs(lhs - rhs.value_) < rhs.epsilon_ + rhs.delta_ * std::abs(rhs.value_);
}

inline bool operator==(const AbsApprox& lhs, double rhs) { return operator==(rhs, lhs); }

inline bool operator!=(double lhs, const AbsApprox& rhs) { return !operator==(lhs, rhs); }

inline bool operator!=(const AbsApprox& lhs, double rhs) { return !operator==(rhs, lhs); }

inline bool operator<=(double lhs, const AbsApprox& rhs) { return lhs < rhs.value_ || lhs == rhs; }

inline bool operator<=(const AbsApprox& lhs, double rhs) { return lhs.value_ < rhs || lhs == rhs; }

inline bool operator>=(double lhs, const AbsApprox& rhs) { return lhs > rhs.value_ || lhs == rhs; }

inline bool operator>=(const AbsApprox& lhs, double rhs) { return lhs.value_ > rhs || lhs == rhs; }

inline bool operator<(double lhs, const AbsApprox& rhs) { return lhs < rhs.value_ && lhs != rhs; }

inline bool operator<(const AbsApprox& lhs, double rhs) { return lhs.value_ < rhs && lhs != rhs; }

inline bool operator>(double lhs, const AbsApprox& rhs) { return lhs > rhs.value_ && lhs != rhs; }

inline bool operator>(const AbsApprox& lhs, double rhs) { return lhs.value_ > rhs && lhs != rhs; }

inline doctest::String toString(const AbsApprox& in) {
  // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
  return doctest::String("AbsApprox( ") + doctest::toString(in.value_) + " )";
}

// test file path derelativizer.
// Used to handle different execution directories for the tests
inline std::string tpath(const std::string& subpath) {
  const auto base = seissol::filesystem::path(TestSetup::Path);
  const auto addend = seissol::filesystem::path(subpath);
  return std::string(base / addend);
}

} // namespace seissol::unit_test

// we add a printer for a vector; as the stdlib doesn't provide one as of C++17

// NOLINTBEGIN (cert-dcl58-cpp)

namespace std {
template <typename T>
ostream& operator<<(ostream& stream, const std::vector<T>& vec) {
  stream << "{ ";
  for (const auto& item : vec) {
    stream << item << ", ";
  }
  stream << "}";
  return stream;
}

template <typename T, std::size_t N>
ostream& operator<<(ostream& stream, const std::array<T, N>& vec) {
  stream << "{ ";
  for (const auto& item : vec) {
    stream << item << ", ";
  }
  stream << "}";
  return stream;
}
} // namespace std

// NOLINTEND

#endif // SEISSOL_TESTS_TESTHELPER_H_
