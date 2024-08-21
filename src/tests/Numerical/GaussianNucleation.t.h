#include "Numerical/GaussianNucleationFunction.h"
#include "doctest.h"
#include <cassert>

namespace seissol::unit_test {
/**
 * Reference implementation for the Gaussian Nucleation function
 */
real gaussianNucleation(real time, real dt, real tau) {
  const auto step = [&tau](real t) {
    if (t < 0) {
      return static_cast<real>(0.0);
    } else if (t > tau) {
      return static_cast<real>(1.0);
    } else {
      return std::exp((t - tau) * (t - tau) / (t * (t - 2 * tau)));
    }
  };
  return step(time) - step(time - dt);
}

TEST_CASE("Gaussian Nucleation Function") {
  constexpr real dt = 0.01;
  constexpr real epsilon = 1e-4;
  for (const real effectiveRiseTime : {0.8, 1.0}) {
    for (int i = -10; i < 111; i++) {
      const real stfEvaluated =
          seissol::gaussianNucleationFunction::smoothStepIncrement(i * dt, dt, effectiveRiseTime);
      const real referenceEvaluated = gaussianNucleation(i * dt, dt, effectiveRiseTime);
      REQUIRE(stfEvaluated == AbsApprox(referenceEvaluated).epsilon(epsilon));
    }
  }
}

} // namespace seissol::unit_test
