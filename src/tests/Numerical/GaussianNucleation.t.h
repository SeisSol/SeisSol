// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/Precision.h"
#include "Numerical/GaussianNucleationFunction.h"
#include "doctest.h"
#include "tests/TestHelper.h"
#include <cassert>

namespace seissol::unit_test {
/**
 * Reference implementation for the Gaussian Nucleation function
 */
template <typename T>
inline T gaussianNucleation(T time, T dt, T tau) {
  const auto step = [&tau](T t) {
    if (t < 0) {
      return static_cast<T>(0.0);
    } else if (t > tau) {
      return static_cast<T>(1.0);
    } else {
      return std::exp((t - tau) * (t - tau) / (t * (t - 2 * tau)));
    }
  };
  return step(time) - step(time - dt);
}

TEST_CASE_TEMPLATE("Gaussian Nucleation Function", RealT, float, double) {
  constexpr RealT Dt = 0.01;
  constexpr RealT Epsilon = 1e-4;
  for (const RealT effectiveRiseTime : {0.8, 1.0}) {
    for (int i = -10; i < 111; i++) {
      const auto stfEvaluated = seissol::gaussianNucleationFunction::smoothStepIncrement<RealT>(
          i * Dt, Dt, effectiveRiseTime);
      const auto referenceEvaluated = gaussianNucleation<RealT>(i * Dt, Dt, effectiveRiseTime);
      REQUIRE(stfEvaluated == AbsApprox(referenceEvaluated).epsilon(Epsilon));
    }
  }
}

} // namespace seissol::unit_test
