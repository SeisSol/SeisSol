// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_PHIFUNCTIONS_H_
#define SEISSOL_SRC_NUMERICAL_PHIFUNCTIONS_H_

#include "Common/Marker.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace seissol::functions {

/**
 * The family of entire functions
 *
 *   phi_n(x) = sum_{k >= 0} x^k / (n+k)!  =  (exp(x) - sum_{k < n} x^k / k!) / x^n
 *
 * and their unnormalized counterparts, the remainders of the exponential series
 *
 *   rem_n(x) = x^n phi_n(x) = exp(x) - sum_{k < n} x^k / k!  =  sum_{k >= n} x^k / k! .
 *
 * The first few remainders are exp(x), exp(x) - 1, exp(x) - 1 - x, exp(x) - 1 - x - x^2/2.
 * phi_n is the quantity that actually occurs in applications, because the division by x^n
 * is usually undone by the surrounding formula; phi_n stays finite at x = 0 (phi_n(0) = 1/n!)
 * while the remainder loses all significant digits there.
 *
 * Evaluating rem_n by its defining difference cancels catastrophically for |x| < n: the
 * result is of size |x|^n/n! while the operands are of size max(1, |x|^{n-1}/(n-1)!).
 * The implementation therefore switches between two forward-stable regimes:
 *
 *   |x| <= n+1: the Taylor series of phi_n. Its terms are non-increasing from the very
 *               first one, so the sum is at least half the largest term and no digits are
 *               lost - for either sign of x. The series is truncated at a compile-time
 *               length chosen for the working precision.
 *   |x| >  n+1: the defining difference. For x > 0 the result is a fixed fraction of
 *               exp(x), for x < 0 the truncated polynomial is dominated by its last term;
 *               either way the operands are no larger than the result.
 *
 * The switching point n+1 is exactly where the two regimes meet, so both branches are
 * accurate to a few ulp everywhere.
 */
namespace phifunctions {

template <typename T>
SEISSOL_HOSTDEVICE constexpr T factorial(std::size_t n) {
  T value{1};
  for (std::size_t i = 2; i <= n; ++i) {
    value *= static_cast<T>(i);
  }
  return value;
}

template <typename T>
SEISSOL_HOSTDEVICE constexpr T integerPower(T x, std::size_t n) {
  T value{1};
  for (std::size_t i = 0; i < n; ++i) {
    value *= x;
  }
  return value;
}

/**
 * Number of series terms after the leading one such that the truncation error stays below
 * the relative tolerance for all |x| <= bound. The terms are measured relative to the
 * leading term 1/n!, which is a lower bound (up to a factor of two) for the sum itself.
 */
constexpr std::size_t seriesLength(std::size_t n, double bound, double tolerance) {
  double term = 1.0;
  std::size_t count = 0;
  while (term > tolerance && count < 1024) {
    ++count;
    term *= bound / static_cast<double>(n + count);
  }
  return count + 1;
}

/// Reciprocals 1/(N+k) of the Horner coefficients, rounded once at compile time.
template <std::size_t N, std::size_t Terms, typename T>
SEISSOL_HOSTDEVICE constexpr std::array<T, Terms + 1> seriesCoefficients() {
  std::array<T, Terms + 1> coefficients{};
  for (std::size_t k = 1; k <= Terms; ++k) {
    coefficients[k] = T{1} / static_cast<T>(N + k);
  }
  return coefficients;
}

/**
 * Horner evaluation of sum_{k = 0}^{Terms} x^k / (N+k)! . Multiplying by the tabulated
 * reciprocals instead of dividing costs at most one further ulp in the sum and turns the
 * dependency chain into Terms fused multiply-adds.
 */
template <std::size_t N, std::size_t Terms, typename T>
SEISSOL_HOSTDEVICE constexpr T phiSeries(T x) {
  constexpr auto Coefficients = seriesCoefficients<N, Terms, T>();
  T value{1};
  for (std::size_t k = Terms; k >= 1; --k) {
    value = T{1} + value * (x * Coefficients[k]);
  }
  return value / factorial<T>(N);
}

/// Horner evaluation of sum_{k < N} x^k / k! .
template <std::size_t N, typename T>
SEISSOL_HOSTDEVICE constexpr T truncatedExponential(T x) {
  if constexpr (N == 0) {
    return T{0};
  } else {
    T value{1};
    for (std::size_t k = N - 1; k >= 1; --k) {
      value = T{1} + value * x / static_cast<T>(k);
    }
    return value;
  }
}

template <std::size_t N, typename T>
constexpr std::size_t DefaultSeriesLength = seriesLength(
    N, static_cast<double>(N + 1), static_cast<double>(std::numeric_limits<T>::epsilon()));

} // namespace phifunctions

/**
 * phi_N(x) = sum_{k >= 0} x^k / (N+k)! , accurate to a few ulp for all finite x.
 *
 * phi_0 = exp(x), phi_1 = expm1(x)/x, phi_2 = (exp(x) - 1 - x)/x^2, ...
 */
template <std::size_t N, typename T = double>
SEISSOL_HOSTDEVICE inline T phi(T x) {
  static_assert(N <= 20, "phi is limited to orders whose factorial is representable");
  if constexpr (N == 0) {
    return std::exp(x);
  } else {
    constexpr T Bound = static_cast<T>(N + 1);
    if (std::abs(x) <= Bound) {
      return phifunctions::phiSeries<N, phifunctions::DefaultSeriesLength<N, T>>(x);
    }
    return (std::exp(x) - phifunctions::truncatedExponential<N, T>(x)) /
           phifunctions::integerPower(x, N);
  }
}

/**
 * phi_N(x) by the Taylor series truncated after Terms+1 terms, with no range check and no
 * exponential. Meant for kernels whose argument range is known to be small, where the
 * generic entry point pays for arguments that never occur and whose branch blocks
 * vectorization. A suitable length follows from a measured bound on |x|:
 *
 *   constexpr auto Terms = phifunctions::seriesLength(N, bound, tolerance);
 *
 * The result is meaningless for |x| beyond that bound, so the bound belongs into an
 * assertion at the call site.
 */
template <std::size_t N, std::size_t Terms, typename T = double>
SEISSOL_HOSTDEVICE constexpr T phiTruncated(T x) {
  return phifunctions::phiSeries<N, Terms, T>(x);
}

/**
 * rem_N(x) = exp(x) - sum_{k < N} x^k / k! , the remainder of the exponential series after
 * N terms, accurate to a few ulp for all finite x.
 *
 * rem_0 = exp(x), rem_1 = expm1(x), rem_2 = exp(x) - 1 - x, rem_3 = exp(x) - 1 - x - x^2/2.
 *
 * For N = 1 this forwards to std::expm1, which is a touch more accurate than the generic
 * path; the two agree to within one ulp.
 */
template <std::size_t N, typename T = double>
SEISSOL_HOSTDEVICE inline T expRemainder(T x) {
  if constexpr (N == 0) {
    return std::exp(x);
  } else if constexpr (N == 1) {
    return std::expm1(x);
  } else {
    return phifunctions::integerPower(x, N) * phi<N, T>(x);
  }
}

namespace phifunctions {
template <typename T, std::size_t... Ns>
SEISSOL_HOSTDEVICE inline std::array<T, sizeof...(Ns)> phiEach(T x,
                                                               std::index_sequence<Ns...> /*idx*/) {
  return {phi<Ns, T>(x)...};
}
} // namespace phifunctions

/**
 * phi_0(x), ..., phi_N(x) in one go, as needed by exponential integrators.
 *
 * For |x| <= 1 a single series evaluation plus the downward recurrence
 * phi_n(x) = x phi_{n+1}(x) + 1/n! suffices; the correction term is then small against
 * 1/n!, so the recurrence is forward-stable. Outside that range every index is evaluated
 * on its own, since the recurrence would subtract nearly equal quantities for small n.
 */
template <std::size_t N, typename T = double>
SEISSOL_HOSTDEVICE inline std::array<T, N + 1> phiUpTo(T x) {
  if (std::abs(x) > T{1}) {
    return phifunctions::phiEach<T>(x, std::make_index_sequence<N + 1>{});
  }
  std::array<T, N + 1> values{};
  values[N] = phi<N, T>(x);
  for (std::size_t n = N; n >= 1; --n) {
    values[n - 1] = x * values[n] + T{1} / phifunctions::factorial<T>(n - 1);
  }
  return values;
}

} // namespace seissol::functions

#endif // SEISSOL_SRC_NUMERICAL_PHIFUNCTIONS_H_
