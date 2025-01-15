// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Numerical/GaussianNucleationFunction.h"
#include "doctest.h"
#include "tests/TestHelper.h"
#include <Kernels/Precision.h>
#include <cassert>

namespace seissol::unit_test {
/**
 * Reference implementation for the Gaussian Nucleation function
 */
inline real gaussianNucleation(real time, real dt, real tau) {
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
  constexpr real Dt = 0.01;
  constexpr real Epsilon = 1e-4;
  for (const real effectiveRiseTime : {0.8, 1.0}) {
    for (int i = -10; i < 111; i++) {
      const real stfEvaluated =
          seissol::gaussianNucleationFunction::smoothStepIncrement(i * Dt, Dt, effectiveRiseTime);
      const real referenceEvaluated = gaussianNucleation(i * Dt, Dt, effectiveRiseTime);
      REQUIRE(stfEvaluated == AbsApprox(referenceEvaluated).epsilon(Epsilon));
    }
  }
}

} // namespace seissol::unit_test
