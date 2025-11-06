// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "doctest.h"

#include "Numerical/RegularizedYoffe.h"

#include <cassert>
#include <cmath>
#include <initializer_list>
#include <limits>

namespace seissol::unit_test {
/**
 * Yoffe function, see Tinti et al. 2005: eq 1
 */
template <typename T>
inline T yoffe(T time, T riseTime) {
  constexpr auto Epsilon = 10 * std::numeric_limits<T>::epsilon();
  if (time < Epsilon || time > riseTime - Epsilon) {
    return 0.0;
  } else {
    return 2 / (M_PI * riseTime) * std::sqrt((riseTime - time) / time);
  }
}

/**
 * Triangle function, see Tinti et al. 2005: eq 2
 */
template <typename T>
inline T triangle(T time, T halfDuration) {
  assert(halfDuration > 0);
  const auto halfDurationSquared = halfDuration * halfDuration;
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
template <typename T>
inline T regularizedYoffe(T time, T tauS, T tauR) {
  if (time < -2 * tauS || time > tauR + 2 * tauS) {
    return 0.0;
  } else {
    // use trapezoidal rule to compute integral
    T yoffeVal = 0.0;
    const size_t numberOfTrapezoids = 1e6;
    // we integrate from -2*tauS until tauR + 2*tauS
    const T h = (tauR + 4 * tauS) / numberOfTrapezoids;
    // integrate yoffe(x) * triangle(time - x) dx
    auto integrand = [&time, &tauS, &tauR](T x) {
      return yoffe(x, tauR) * triangle(time - x, tauS);
    };
    // integrand is 0 at left and right boundary, only accumulate interior values
    for (size_t i = 1; i < numberOfTrapezoids; ++i) {
      yoffeVal += integrand(i * h);
    }
    return yoffeVal * h;
  }
}

TEST_CASE_TEMPLATE("Regularized Yoffe Function", RealT, float, double) {
  constexpr RealT Dt = 0.01;
  // rather coarse epsilon, since the quadrature to compute the reference is also rather coarse
  constexpr RealT Epsilon = 1e-2;
  for (const RealT accTime : {0.2, 0.3}) {
    for (const RealT effectiveRiseTime : {0.9, 1.1}) {
      const RealT tauS = accTime / 1.27;
      const RealT tauR = effectiveRiseTime - 2 * tauS;

      for (int i = -10; i < 111; i++) {
        const auto stfEvaluated =
            seissol::regularizedYoffe::regularizedYoffe<RealT>(i * Dt, tauS, tauR);
        const auto referenceEvaluated = regularizedYoffe<RealT>(i * Dt, tauS, tauR);
        REQUIRE(stfEvaluated == AbsApprox(referenceEvaluated).epsilon(Epsilon));
      }
    }
  }
}

} // namespace seissol::unit_test
