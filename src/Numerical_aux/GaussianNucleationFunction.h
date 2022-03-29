#ifndef SEISSOL_GAUSSIANNUCELATIONFUNCTION_H
#define SEISSOL_GAUSSIANNUCELATIONFUNCTION_H

namespace seissol::gaussianNucleationFunction {
/**
 * Implementation of the gaussian nucleation function, which is widely used in the TPV benchmarks.
 */

/**
 * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
 */
inline real smoothStep(real currentTime, real t0) {
  if (currentTime <= 0) {
    return 0.0;
  } else if (currentTime < t0) {
    real tau = currentTime - t0;
    return std::exp(tau * tau /
                    (currentTime * (currentTime - 2.0 * t0)));
  } else {
    return 1.0;
  }
}

/**
 * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
 */
inline real smoothStepIncrement(real currentTime, real dt, real t0) {
  real gNuc = smoothStep(currentTime, t0);
  real prevTime = currentTime - dt;
  gNuc = gNuc - smoothStep(prevTime, t0);
  return gNuc;
}

} // namespace seissol::gaussianNucleationFunction

#endif // SEISSOL_GAUSSIANNUCELATIONFUNCTION_H
