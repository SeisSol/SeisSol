// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_TESTHELPER_H_
#define SEISSOL_TESTS_TESTHELPER_H_

#include "doctest.h"

#include <cmath>
#include <limits>
#include <ostream>
#include <vector>

namespace seissol::unit_test {

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
  double m_epsilon{std::numeric_limits<double>::epsilon()};
  double m_delta{0.0};
  double m_value{0.0};
};

AbsApprox::AbsApprox(double value) : m_value(value) {}

AbsApprox AbsApprox::operator()(double newValue) const {
  AbsApprox approx(newValue);
  approx.epsilon(m_epsilon);
  approx.delta(m_delta);
  return approx;
}

AbsApprox& AbsApprox::epsilon(double newEpsilon) {
  m_epsilon = newEpsilon;
  return *this;
}

AbsApprox& AbsApprox::delta(double newDelta) {
  m_delta = newDelta;
  return *this;
}

inline bool operator==(double lhs, const AbsApprox& rhs) {
  return std::abs(lhs - rhs.m_value) < rhs.m_epsilon + rhs.m_delta * std::abs(rhs.m_value);
}

inline bool operator==(const AbsApprox& lhs, double rhs) { return operator==(rhs, lhs); }

inline bool operator!=(double lhs, const AbsApprox& rhs) { return !operator==(lhs, rhs); }

inline bool operator!=(const AbsApprox& lhs, double rhs) { return !operator==(rhs, lhs); }

inline bool operator<=(double lhs, const AbsApprox& rhs) { return lhs < rhs.m_value || lhs == rhs; }

inline bool operator<=(const AbsApprox& lhs, double rhs) { return lhs.m_value < rhs || lhs == rhs; }

inline bool operator>=(double lhs, const AbsApprox& rhs) { return lhs > rhs.m_value || lhs == rhs; }

inline bool operator>=(const AbsApprox& lhs, double rhs) { return lhs.m_value > rhs || lhs == rhs; }

inline bool operator<(double lhs, const AbsApprox& rhs) { return lhs < rhs.m_value && lhs != rhs; }

inline bool operator<(const AbsApprox& lhs, double rhs) { return lhs.m_value < rhs && lhs != rhs; }

inline bool operator>(double lhs, const AbsApprox& rhs) { return lhs > rhs.m_value && lhs != rhs; }

inline bool operator>(const AbsApprox& lhs, double rhs) { return lhs.m_value > rhs && lhs != rhs; }

inline doctest::String toString(const AbsApprox& in) {
  // NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks)
  return doctest::String("AbsApprox( ") + doctest::toString(in.m_value) + " )";
}

} // namespace seissol::unit_test

// we add a printer for a vector; as the stdlib doesn't provide one as of C++17

// NOLINTBEGIN (cert-dcl58-cpp)

namespace std {
template <typename T>
ostream& operator<<(ostream& stream, const std::vector<T>& vec) {
  stream << "{";
  for (const auto& item : vec) {
    stream << item << ",";
  }
  stream << "}";
  return stream;
}

template <typename T, std::size_t N>
ostream& operator<<(ostream& stream, const std::array<T, N>& vec) {
  stream << "{";
  for (const auto& item : vec) {
    stream << item << ",";
  }
  stream << "}";
  return stream;
}
} // namespace std

// NOLINTEND

#endif // SEISSOL_TESTS_TESTHELPER_H_
