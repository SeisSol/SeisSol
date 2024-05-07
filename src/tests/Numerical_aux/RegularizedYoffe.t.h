#include "doctest.h"
#include "Numerical_aux/RegularizedYoffe.h"
#include <cassert>

namespace seissol::unit_test {
/**
 * Yoffe function, see Tinti et al. 2005: eq 1
 */
real yoffe(real time, real riseTime) {
  constexpr real epsilon = 10 * std::numeric_limits<real>::epsilon();
  if (time < epsilon || time > riseTime - epsilon) {
    return 0.0;
  } else {
    return 2 / (M_PI * riseTime) * std::sqrt((riseTime - time) / time);
  }
}

/**
 * Triangle function, see Tinti et al. 2005: eq 2
 */
real triangle(real time, real halfDuration) {
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
real regularizedYoffe(real time, real tauS, real tauR) {
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
  constexpr real dt = 0.01;
  // rather coarse epsilon, since the quadrature to compute the reference is also rather coarse
  constexpr real epsilon = 1e-2;
  for (real const accTime : {0.2, 0.3}) {
    for (real const effectiveRiseTime : {0.9, 1.1}) {
      real tauS = accTime / 1.27;
      real tauR = effectiveRiseTime - 2 * tauS;

      for (int i = -10; i < 111; i++) {
        real stfEvaluated = seissol::regularizedYoffe::regularizedYoffe(i * dt, tauS, tauR);
        real referenceEvaluated = regularizedYoffe(i * dt, tauS, tauR);
        REQUIRE(stfEvaluated == AbsApprox(referenceEvaluated).epsilon(epsilon));
      }
    }
  }
}

} // namespace seissol::unit_test