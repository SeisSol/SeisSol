// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATECOMMON_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATECOMMON_H_

#include "Common/Marker.h"
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

struct Settings {
  /**
   * Parameters of the optimisation loops
   * absolute tolerance on the function to be optimized
   * This value is quite arbitrary (a bit bigger as the expected numerical error) and may not be
   * the most adapted Number of iteration in the loops
   */

  const uint32_t maxNumberSlipRateUpdates{60};
  const uint32_t numberStateVariableUpdates{2};
  const double newtonTolerance{1e-8};
};

/**
  Computes asinh(x * exp(c)). Reason is: exp(c) can grow really large (too large for float);
  but actually asinh(exp(c)) \approx c for large c.

  Hence, we compute instead
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
  if (expLog > 0) {
    return expLog + std::log(x + std::sqrt(x * x + exp));
  } else {
    const auto v = exp * x;
    return std::asinh(v);
  }
}

/**
  Helper function to arsinhexp. Since for asinh(x * exp(c)),
  we can assume c to be constant, we can pre-compute exp(c) or exp(-2c).
 */
#pragma omp declare simd
template <typename T>
SEISSOL_HOSTDEVICE constexpr T computeCExp(T cExpLog) {
  T cExp{};
  if (cExpLog > 0) {
    cExp = std::exp(-2 * cExpLog);
  } else {
    cExp = std::exp(cExpLog);
  }
  return cExp;
}

/**
  Derivative to arsinhexp.
 */
#pragma omp declare simd
template <typename T>
SEISSOL_HOSTDEVICE constexpr T arsinhexpDerivative(T x, T expLog, T exp) {
  if (expLog > 0) {
    return 1 / std::sqrt(x * x + exp);
  } else {
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
  const T sign = c > 0 ? 1 : -1;
  const T absC = std::abs(c);
  return absC + std::log(x / 2 * -sign * std::expm1(-2 * absC));
}

} // namespace seissol::dr::friction_law::rs

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATECOMMON_H_
