// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Numerical/RegularizedYoffe.h"
#include "doctest.h"
#include <cassert>

namespace seissol::unit_test {
/**
 * Yoffe function, see Tinti et al. 2005: eq 1
 */
inline real yoffe(real time, real riseTime) {
  constexpr real Epsilon = 10 * std::numeric_limits<real>::epsilon();
  if (time < Epsilon || time > riseTime - Epsilon) {
    return 0.0;
  } else {
    return 2 / (M_PI * riseTime) * std::sqrt((riseTime - time) / time);
  }
}

/**
 * Triangle function, see Tinti et al. 2005: eq 2
 */
inline real triangle(real time, real halfDuration) {
  assert(halfDuration > 0);
  const real halfDurationSquared = halfDuration * halfDuration;
  if (time > 0 && time < halfDuration) {
    return time / halfDurationSquared;
  } else if (time > halfDuration && time < 2 * halfDuration) {
    return (2 * halfDuration - time) / halfDurationSquared;
  } else {
    return 0;
  }
}

/**
 * Reference implementation for the regularized Yoffe function by explicitly computing the
 * convolution, see Tinti et al. 2005: eq 3
 */
inline real regularizedYoffe(real time, real tauS, real tauR) {
  if (time < -2 * tauS || time > tauR + 2 * tauS) {
    return 0.0;
  } else {
    // use trapezoidal rule to compute integral
    real yoffeVal = 0.0;
    const size_t numberOfTrapezoids = 1e6;
    // we integrate from -2*tauS until tauR + 2*tauS
    const real h = (tauR + 4 * tauS) / numberOfTrapezoids;
    // integrate yoffe(x) * triangle(time - x) dx
    auto integrand = [&time, &tauS, &tauR](real x) {
      return yoffe(x, tauR) * triangle(time - x, tauS);
    };
    // integrand is 0 at left and right boundary, only accumulate interior values
    for (size_t i = 1; i < numberOfTrapezoids; i++) {
      yoffeVal += integrand(i * h);
    }
    return yoffeVal * h;
  }
}

TEST_CASE("Regularized Yoffe Function") {
  constexpr real Dt = 0.01;
  // rather coarse epsilon, since the quadrature to compute the reference is also rather coarse
  constexpr real Epsilon = 1e-2;
  for (const real accTime : {0.2, 0.3}) {
    for (const real effectiveRiseTime : {0.9, 1.1}) {
      const real tauS = accTime / 1.27;
      const real tauR = effectiveRiseTime - 2 * tauS;

      for (int i = -10; i < 111; i++) {
        real stfEvaluated = seissol::regularizedYoffe::regularizedYoffe(i * Dt, tauS, tauR);
        const real referenceEvaluated = regularizedYoffe(i * Dt, tauS, tauR);
        REQUIRE(stfEvaluated == AbsApprox(referenceEvaluated).epsilon(Epsilon));
      }
    }
  }
}

} // namespace seissol::unit_test
