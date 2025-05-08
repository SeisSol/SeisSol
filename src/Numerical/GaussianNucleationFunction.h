// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_GAUSSIANNUCLEATIONFUNCTION_H_
#define SEISSOL_SRC_NUMERICAL_GAUSSIANNUCLEATIONFUNCTION_H_

#include "Numerical/Functions.h"

namespace seissol::gaussianNucleationFunction {
/**
 * Implementation of the gaussian nucleation function, which is widely used in the TPV benchmarks.
 */

/**
 * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
 */
template <typename T, typename MathFunctions>
SEISSOL_HOSTDEVICE inline T smoothStep(T currentTime, T t0) {
  if (currentTime <= 0) {
    return 0.0;
  } else if (currentTime < t0) {
    const T tau = currentTime - t0;
    return MathFunctions::exp(tau * tau / (currentTime * (currentTime - 2.0 * t0)));
  } else {
    return 1.0;
  }
}

/**
 * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
 */
template <typename T, typename MathFunctions = seissol::functions::HostStdFunctions>
SEISSOL_HOSTDEVICE inline T smoothStepIncrement(T currentTime, T dt, T t0) {
  return smoothStep<T, MathFunctions>(currentTime, t0) -
         smoothStep<T, MathFunctions>(currentTime - dt, t0);
}

} // namespace seissol::gaussianNucleationFunction

#endif // SEISSOL_SRC_NUMERICAL_GAUSSIANNUCLEATIONFUNCTION_H_
