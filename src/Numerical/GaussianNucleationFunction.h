// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_GAUSSIANNUCLEATIONFUNCTION_H_
#define SEISSOL_SRC_NUMERICAL_GAUSSIANNUCLEATIONFUNCTION_H_

#include "Numerical/Functions.h"

#include <cmath>

namespace seissol::gaussianNucleationFunction {
/**
 * Implementation of the gaussian nucleation function, which is widely used in the TPV benchmarks.
 */

/**
 * Exponent of the nucleation function on the ramp, @f$ (t - t_0)^2 / (t (t - 2 t_0)) @f$.
 *
 * Of the algebraically equal ways of writing it, this is the one that stays accurate near
 * @f$ t = t_0 @f$; the form @f$ 1 + t_0^2 / (t (t - 2 t_0)) @f$ adds two nearly opposite
 * numbers there. It is the other way round for the difference of two exponents, see below.
 */
template <typename T>
SEISSOL_HOSTDEVICE inline T smoothStepExponent(T currentTime, T t0) {
  const T tau = currentTime - t0;
  return tau * tau / (currentTime * (currentTime - static_cast<T>(2) * t0));
}

/**
 * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
 */
template <typename T>
SEISSOL_HOSTDEVICE inline T smoothStep(T currentTime, T t0) {
  if (currentTime <= 0) {
    return static_cast<T>(0);
  } else if (currentTime < t0) {
    return std::exp(smoothStepExponent(currentTime, t0));
  } else {
    return static_cast<T>(1);
  }
}

/**
 * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
 *
 * The two values agree to within dt / t0 over most of the ramp, so subtracting them gives up
 * that many digits, and the caller sums the results over the whole nucleation. Both are
 * exponentials, so the difference factors into @f$ f(t-dt) (e^{g(t) - g(t-dt)} - 1) @f$ with
 * an expm1 as the second factor. The difference of the exponents has a closed form that
 * subtracts nothing nearly equal either: from @f$ g(t) = 1 + t_0^2 / (t (t - 2 t_0)) @f$,
 * @f[ g(t) - g(t - dt) =
 *     \frac{t_0^2 dt (2 t_0 + dt - 2 t)}{t (t - 2 t_0) (t - dt) (t - dt - 2 t_0)} @f]
 *
 * Where the two values are an e-fold or more apart, the factored form would multiply an
 * underflowed exponential by an overflowed expm1, and the plain difference is the one that
 * is well conditioned.
 */
template <typename T>
SEISSOL_HOSTDEVICE inline T smoothStepIncrement(T currentTime, T dt, T t0) {
  const T previousTime = currentTime - dt;
  if (currentTime <= 0 || previousTime >= t0) {
    return static_cast<T>(0);
  }
  if (previousTime <= 0) {
    return smoothStep<T>(currentTime, t0);
  }
  const T previousExponent = smoothStepExponent(previousTime, t0);
  if (currentTime >= t0) {
    return -std::expm1(previousExponent);
  }
  const T currentScale = currentTime * (currentTime - static_cast<T>(2) * t0);
  const T previousScale = previousTime * (previousTime - static_cast<T>(2) * t0);
  // the step the two exponents are actually a distance apart, which is not dt once
  // currentTime - dt has been rounded, and an error here would show up divided by dt
  const T step = currentTime - previousTime;
  // step + 2 (t0 - t) rather than 2 t0 + step - 2 t: near the end of the ramp the latter
  // rounds a sum of size t0 and then subtracts it down to a result of size step
  const T deltaExponent = t0 * t0 * step * (step + static_cast<T>(2) * (t0 - currentTime)) /
                          (currentScale * previousScale);
  if (deltaExponent >= static_cast<T>(1)) {
    return std::exp(smoothStepExponent(currentTime, t0)) - std::exp(previousExponent);
  }
  return std::exp(previousExponent) * std::expm1(deltaExponent);
}

} // namespace seissol::gaussianNucleationFunction

#endif // SEISSOL_SRC_NUMERICAL_GAUSSIANNUCLEATIONFUNCTION_H_
