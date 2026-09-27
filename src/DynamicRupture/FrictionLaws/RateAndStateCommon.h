// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATECOMMON_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATECOMMON_H_

#include "Common/Marker.h"
#include "DynamicRupture/FrictionLaws/Dual.h"
#include "Kernels/Precision.h"

#include <cstdint>
#include <type_traits>

namespace seissol::dr::friction_law::rs {
// If the SR is too close to zero, we will have problems (NaN)
// as a consequence, the SR is affected the AlmostZero value when too small
// For double precision 1e-45 is a chosen by trial and error. For single precision, this value is
// too small, so we use 1e-35
constexpr real almostZero() {
  if constexpr (std::is_same<real, double>()) {
    return 1e-45;
  } else if constexpr (std::is_same<real, float>()) {
    return 1e-35;
  } else {
    return std::numeric_limits<real>::min();
  }
}

/**
  Computes asinh(x * exp(c)). Reason is: exp(c) can grow really large (too large for float);
  but actually asinh(exp(c)) \approx c for large c.

  Hence, we compute instead (x > 0)
  asinh(x * exp(c))
  = asinh((x * exp(c)) + sqrt((x * exp(c))**2 + 1))
  = asinh(exp(c) * (x + sqrt(x**2 + exp(-2c))))
  = c + asinh(x + sqrt(x**2 + exp(-2c))).

  Here, exp(-2c) is small.

  If c < 0, we can process as normal.
 */
#pragma omp declare simd
template <typename T>
SEISSOL_HOSTDEVICE constexpr T arsinhexp(T x, T expLog, T exp) {
  // unqualified so that a dual number picks up the overloads next to its own definition, while a
  // plain scalar keeps the standard ones
  using std::abs;
  using std::asinh;
  using std::log;
  using std::sqrt;
  using Scalar = decltype(valueOf(T{}));

  // Switch is empirically chosen; to prevent issues with
  // or replacement formula not being accurate enough if x * exp(c) is small
  constexpr Scalar Switch = 10;
  constexpr Scalar Threshold = 50;
  constexpr Scalar Log2 = 0.69314718055994530943;
  int xexp{};
  (void)std::frexp(valueOf(x), &xexp);

  // make sure to invert the constant we'd use otherwise (if the exponent is too big/small)

  // the branch selects a formula; the selected formula is what carries the derivative
  // use the new code path only if we really need to
  if (valueOf(expLog) + std::max(xexp, 0) * Log2 > Switch || valueOf(expLog) >= Threshold) {
    if (valueOf(expLog) <= 0) {
      exp = T(1) / exp;
    }
    const T xa = abs(x);
    const T xs = valueOf(x) >= 0 ? T(1) : T(-1);
    return xs * (expLog + log(xa + sqrt(xa * xa + exp * exp)));
  } else {
    if (valueOf(expLog) > 0) {
      exp = T(1) / exp;
    }
    const auto v = exp * x;
    return asinh(v);
  }
}

/**
  Helper function to arsinhexp. Since for asinh(x * exp(c)),
  we can assume c to be constant, we can pre-compute exp(c) or exp(-2c).
 */
#pragma omp declare simd
template <typename T>
SEISSOL_HOSTDEVICE constexpr T computeCExp(T cExpLog) {
  using std::exp;
  T cExp{};
  if (valueOf(cExpLog) > 0) {
    cExp = exp(-cExpLog);
  } else {
    cExp = exp(cExpLog);
  }
  return cExp;
}

/**
  Derivative to arsinhexp.
 */
#pragma omp declare simd
template <typename T>
SEISSOL_HOSTDEVICE constexpr T arsinhexpDerivative(T x, T expLog, T exp) {
  constexpr T Switch = 10;
  constexpr T Threshold = 50;
  constexpr T Log2 = 0.69314718055994530943;
  int xexp{};
  (void)std::frexp(x, &xexp);

  // make sure to invert the constant we'd use otherwise (if the exponent is too big/small)

  if (expLog + std::max(xexp, 0) * Log2 > Switch || expLog >= Threshold) {
    if (expLog <= 0) {
      exp = 1 / exp;
    }
    return 1 / std::sqrt(x * x + exp * exp);
  } else {
    if (expLog > 0) {
      exp = 1 / exp;
    }
    const auto v = exp * x;
    return exp / std::sqrt(1 + v * v);
  }
}

/**
  Compute log(x * sinh(c)).
  if c > 0, then
  log(x * (e(c) - e(-c)) / 2)
  = log(e(c) * x / 2 * (1 - e(-2c)))
  = c + log(x / 2 * -expm1(-2c))

  if c < 0, then
  log(x * sinh(c)) = -c + log(x / 2 * expm1(2c))

  In total,
  log(x * sinh(c)) = |c| + log(x / 2 * -sign(c) * expm1(-2|c|))
 */
#pragma omp declare simd
template <typename T>
SEISSOL_HOSTDEVICE constexpr T logsinh(T x, T c) {
  using std::abs;
  using std::expm1;
  using std::log;
  const T sign = valueOf(c) >= 0 ? T(1) : T(-1);
  const T absC = abs(c);
  return absC + log(x / T(2) * -sign * expm1(T(-2) * absC));
}

} // namespace seissol::dr::friction_law::rs

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATECOMMON_H_
