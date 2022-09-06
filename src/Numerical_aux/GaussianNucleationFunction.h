#ifndef SEISSOL_GAUSSIANNUCELATIONFUNCTION_H
#define SEISSOL_GAUSSIANNUCELATIONFUNCTION_H

#include "Numerical_aux/Functions.h"

namespace seissol::gaussianNucleationFunction {
/**
 * Implementation of the gaussian nucleation function, which is widely used in the TPV benchmarks.
 */


/**
 * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
 */
template <typename MathFunctions>
inline real smoothStep(real currentTime, real t0) {
  if (currentTime <= 0) {
    return 0.0;
  } else if (currentTime < t0) {
    const real tau = currentTime - t0;
    return MathFunctions::exp(tau * tau /
                    (currentTime * (currentTime - 2.0 * t0)));
  } else {
    return 1.0;
  }
}

/**
 * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
 */
template <typename MathFunctions = seissol::functions::HostStdFunctions>
inline real smoothStepIncrement(real currentTime, real dt, real t0) {
  return smoothStep<MathFunctions>(currentTime, t0) - smoothStep<MathFunctions>(currentTime - dt, t0);
}

} // namespace seissol::gaussianNucleationFunction

#endif // SEISSOL_GAUSSIANNUCELATIONFUNCTION_H
