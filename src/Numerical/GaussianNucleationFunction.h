// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_GAUSSIANNUCLEATIONFUNCTION_H_
#define SEISSOL_SRC_NUMERICAL_GAUSSIANNUCLEATIONFUNCTION_H_

#include "Kernels/Precision.h"
#include "Numerical/Functions.h"

namespace seissol::gaussianNucleationFunction {
/**
 * Implementation of the gaussian nucleation function, which is widely used in the TPV benchmarks.
 */

/**
 * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
 */
#pragma omp declare simd
template <typename MathFunctions>
SEISSOL_HOSTDEVICE inline real smoothStep(real currentTime, real t0) {
  if (currentTime <= 0) {
    return 0.0;
  } else if (currentTime < t0) {
    const real tau = currentTime - t0;
    return MathFunctions::exp(tau * tau / (currentTime * (currentTime - 2.0 * t0)));
  } else {
    return 1.0;
  }
}

// start situation: G(T - dt) = 0; i.e. we have (G(T) - 0 = G(T)) as increment
#pragma omp declare simd
template <typename MathFunctions>
SEISSOL_HOSTDEVICE inline real smoothStepIncrementBegin(real currentTime, real t0) {
  const real tau = currentTime - t0;
  return MathFunctions::exp(tau * tau / (currentTime * (currentTime - 2.0 * t0)));
}

// end situation: G(T) = 1; thus we obtain (1 - G(T - dt)) as increment
#pragma omp declare simd
template <typename MathFunctions>
SEISSOL_HOSTDEVICE inline real smoothStepIncrementEnd(real currentTime, real t0) {
  const real tau = currentTime - t0;
  return -MathFunctions::expm1(tau * tau / (currentTime * (currentTime - 2.0 * t0)));
}

// all other situations: we have a function of the form (exp(a) - exp(b)) as increment
// reformulate (exp(a) - exp(b)) = exp(b) * (exp(a - b) - 1)
#pragma omp declare simd
template <typename MathFunctions>
SEISSOL_HOSTDEVICE inline real smoothStepIncrementMiddle(real currentTime, real dt, real t0) {
  const real previousTime = currentTime - dt;
  const real tau0 = previousTime - t0;
  const real tau1 = currentTime - t0;

  const real valA = tau1 * tau1 / (currentTime * (currentTime - 2.0 * t0));
  const real valB = tau0 * tau0 / (previousTime * (previousTime - 2.0 * t0));

  const real expB = MathFunctions::exp(valB);
  const real expm1AB = MathFunctions::expm1(valA - valB);

  return expB * expm1AB;
}

/**
 * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
 */
#pragma omp declare simd
template <typename MathFunctions = seissol::functions::HostStdFunctions>
SEISSOL_HOSTDEVICE inline real smoothStepIncrement(real currentTime, real dt, real t0) {
  /*
   * We compute
  return smoothStep<MathFunctions>(currentTime, t0) -
         smoothStep<MathFunctions>(currentTime - dt, t0);
  */
  if (currentTime <= 0) {
    return 0.0;
  } else if (currentTime < std::min(t0, dt)) {
    return smoothStepIncrementBegin<MathFunctions>(currentTime, t0);
  } else if (currentTime >= t0 && currentTime < t0 + dt) {
    return smoothStepIncrementEnd<MathFunctions>(currentTime - dt, t0);
  } else if (currentTime < t0) {
    // implicitly by else if 1: currentTime >= dt
    return smoothStepIncrementMiddle<MathFunctions>(currentTime, dt, t0);
  } else {
    return 0.0;
  }
}

} // namespace seissol::gaussianNucleationFunction

#endif // SEISSOL_SRC_NUMERICAL_GAUSSIANNUCLEATIONFUNCTION_H_
