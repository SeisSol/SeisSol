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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
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
  The largest argument whose exponential is comfortably representable. The slack to log(max()) --
  88.7 for float, 709.8 for double -- absorbs the one binade by which the bound on |x| in
  arsinhexp may overshoot.
 */
template <typename T>
SEISSOL_HOSTDEVICE constexpr T logMaxExp() {
  return std::is_same_v<T, float> ? T(87) : T(700);
}

/**
  Precomputes exp(c) for arsinhexp. c does not depend on the slip rate, so for a friction law whose
  state variable stays outside the inversion this runs once per point and time step.

  Returns zero where exp(c) is not representable; arsinhexp then takes its asymptotic branch and
  never reads the value. Zero rather than infinity is deliberate: a masked SIMD loop evaluates both
  branches on every lane, and inf * 0 raises FE_INVALID on a locked point where 0 * 0 stays quiet.
  It also survives the licence -ffast-math grants the compiler to assume that no infinities exist.
 */
#pragma omp declare simd
template <typename T>
SEISSOL_HOSTDEVICE constexpr T computeCExp(T cExpLog) {
  // unqualified so that a dual number picks up the overload next to its own definition, while a
  // plain scalar keeps the standard one
  using std::exp;
  using Scalar = decltype(valueOf(T{}));
  return valueOf(cExpLog) < logMaxExp<Scalar>() ? exp(cExpLog) : T(0);
}

/**
  Computes asinh(x * exp(c)), with c = cExpLog and cExp = exp(c) precomputed by computeCExp. The
  point is that exp(c) alone overflows long before asinh(x * exp(c)) does -- c reaches a few
  hundred for a locked point, while the result stays of the order of c itself.

  frexp bounds log2|x| without evaluating a logarithm: |x| < 2^xexp, hence
  cExp * x < exp(c + xexp * log 2). Clamping the exponent at zero also forces c < logMaxExp, which
  covers |x| < 1, where exp(c) alone is the binding constraint. So the test never admits a product
  that overflows, and wherever the product is representable the plain formula is what runs.

  Where it is not, x * exp(c) lies far beyond 1 / sqrt(eps), and there asinh(z) = log(2z) holds to
  machine precision -- the asymptotic branch therefore needs neither exp nor asinh. It is odd in x,
  like asinh itself, and returns zero at x = 0: the asymptotic form has a logarithmic singularity
  there which the function it stands in for does not.

  The two branches cover everything reachable from a friction law, where x = V / (2 V_0) with V
  clamped from below by almostZero(): the asymptotic branch is then only ever entered at a product
  above 1e8, decades beyond where it becomes exact. Two regions outside that are inaccurate, and a
  caller stepping outside should know which. Where exp(c) overflows while |x| is small enough to
  bring the product back into range -- below 1e-38 in single precision -- the asymptotic branch is
  entered at a product of order one and is simply the wrong formula. Where exp(c) underflows while
  |x| is large enough to lift the product back, mirroring the first, the precomputed factor is zero
  and the plain branch returns zero. Both would need exp(c/2) and two multiplications in place of
  one, which costs a sixth of this function in the folded inversion -- measured -- for accuracy
  outside the domain the friction laws occupy.
 */
#pragma omp declare simd
template <typename T>
SEISSOL_HOSTDEVICE constexpr T arsinhexp(T x, T cExpLog, T cExp) {
  using std::abs;
  using std::asinh;
  using std::log;
  using Scalar = decltype(valueOf(T{}));
  constexpr Scalar Log2 = 0.69314718055994530943;

  int xexp{};
  (void)std::frexp(valueOf(x), &xexp);

  // the branch selects a formula; the selected formula is what carries the derivative
  if (valueOf(cExpLog) + std::max(xexp, 0) * Log2 < logMaxExp<Scalar>()) {
    return asinh(cExp * x);
  }
  if (valueOf(x) == 0) {
    return T(0);
  }
  const T xs = valueOf(x) >= 0 ? T(1) : T(-1);
  return xs * (T(Log2) + cExpLog + log(abs(x)));
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
