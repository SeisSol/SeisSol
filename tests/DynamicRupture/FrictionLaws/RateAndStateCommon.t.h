// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATECOMMON_T_H_
#define SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATECOMMON_T_H_

#include <doctest.h>

#include "DynamicRupture/FrictionLaws/Dual.h"
#include "DynamicRupture/FrictionLaws/RateAndStateCommon.h"
#include "Kernels/Precision.h"
#include "TestHelper.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

namespace seissol::unit_test {

namespace rstest {

namespace rs = seissol::dr::friction_law::rs;
using seissol::dr::friction_law::Dual;

/**
 * Reference for asinh(x exp(c)), in double precision and sharing no branch with the function under
 * test.
 *
 * Where the product is representable, the naive expression is the reference. Beyond that,
 *     asinh(z) = log(z) + log(1 + sqrt(1 + z^-2)),
 * which is exact for every positive z and needs only the logarithm of the product, c + log|x|.
 * The correction term settles at log 2 but is carried rather than assumed, so the reference also
 * says how much of the result the asymptotic form leaves out.
 */
inline double reference(double x, double c) {
  const double product = x * std::exp(c);
  if (std::abs(c) < 700.0 && std::isfinite(product)) {
    return std::asinh(product);
  }
  if (x == 0.0) {
    return 0.0;
  }
  const double logProduct = c + std::log(std::abs(x));
  const double magnitude = logProduct + std::log1p(std::sqrt(1.0 + std::exp(-2.0 * logProduct)));
  return x >= 0.0 ? magnitude : -magnitude;
}

inline double relativeError(double value, double target) {
  if (target == 0.0) {
    return std::abs(value);
  }
  return std::abs((value - target) / target);
}

/// The calling convention of every call site: the caller precomputes exp(c) and passes both.
template <typename T>
T evaluate(T x, T cExpLog) {
  return rs::arsinhexp(x, cExpLog, rs::computeCExp(cExpLog));
}

/// The floor almostZero() puts under the slip rate in a build of precision T. Checked against the
/// function itself for whichever precision this build uses.
template <typename T>
constexpr double slipRateFloor() {
  return std::is_same_v<T, float> ? 1e-35 : 1e-45;
}

/**
 * The slip rates a friction law can present, times the 1 / (2 V_0) the laws fold into x. Below the
 * floor the inversion does not go, and above a hundred metres per second there is no fault.
 */
template <typename T>
std::vector<double> reachableX() {
  const auto scale = 0.5 / 1e-6;
  std::vector<double> values;
  for (double exponent = std::log10(scale * slipRateFloor<T>());
       exponent <= std::log10(scale * 100.0);
       exponent += 0.5) {
    values.push_back(std::pow(10.0, exponent));
  }
  return values;
}

/**
 * Psi / a, which reaches a few hundred at a locked point and turns negative nowhere physical -- the
 * fault initializer warns where it finds a negative state variable.
 *
 * The lower end stops where exp(c) stops being a normal number of T, because from there the
 * precomputed factor is itself a subnormal and carries its own error into everything downstream.
 * That is one of the two bands the header names, and it is checked as such further down.
 */
template <typename T>
std::vector<double> sampleCExpLog() {
  const double floor = std::log(static_cast<double>(std::numeric_limits<T>::min()));
  std::vector<double> values;
  for (double c = floor + 2.5; c <= 400.0; c += 2.5) {
    values.push_back(c);
  }
  return values;
}

} // namespace rstest

// ---------------------------------------------------------------------------
// almostZero
// ---------------------------------------------------------------------------

TEST_CASE("DR RateAndState almostZero" * doctest::test_suite("dynamicrupture")) {
  const auto value = rstest::rs::almostZero();

  CHECK(value > static_cast<real>(0.0));
  CHECK(std::isfinite(value));

  // the floor the inversion clamps its bracket to. A subnormal -- or a value flushed to zero by
  // -ffast-math -- would put the floor back at zero and the logarithm of the slip rate at minus
  // infinity.
  CHECK(std::isnormal(value));
  CHECK(value < static_cast<real>(1e-30));

  if constexpr (std::is_same_v<real, double>) {
    CHECK(value == static_cast<real>(1e-45));
  } else if constexpr (std::is_same_v<real, float>) {
    CHECK(value == static_cast<real>(1e-35));
  }
  // the domain the accuracy checks below sweep starts here
  CHECK(value == static_cast<real>(rstest::slipRateFloor<real>()));
}

// ---------------------------------------------------------------------------
// logMaxExp
// ---------------------------------------------------------------------------

TEST_CASE("DR RateAndState logMaxExp leaves room for one binade" *
          doctest::test_suite("dynamicrupture")) {
  using rstest::rs::logMaxExp;

  SUBCASE("Its exponential is representable") {
    CHECK(std::isfinite(std::exp(logMaxExp<float>())));
    CHECK(std::isfinite(std::exp(logMaxExp<double>())));
  }

  SUBCASE("With a factor of two to spare") {
    // the bound on |x| in arsinhexp comes from frexp and may overshoot the true magnitude by up to
    // one binade, so the product may be twice what the guard admitted
    CHECK(std::isfinite(2.0F * std::exp(logMaxExp<float>())));
    CHECK(std::isfinite(2.0 * std::exp(logMaxExp<double>())));
  }

  SUBCASE("And not more room than it needs") {
    // a threshold far below the range would push the asymptotic branch into arguments where it is
    // not yet exact
    CHECK(logMaxExp<float>() > 80.0F);
    CHECK(logMaxExp<double>() > 690.0);
  }
}

// ---------------------------------------------------------------------------
// computeCExp
// ---------------------------------------------------------------------------

TEST_CASE("DR RateAndState computeCExp" * doctest::test_suite("dynamicrupture")) {
  using rstest::Dual;
  using rstest::rs::computeCExp;
  using rstest::rs::logMaxExp;

  SUBCASE("The exponential itself where it is representable") {
    for (const double c : {-300.0, -120.0, -50.0, -10.0, -1.0, 0.0, 1.0, 10.0, 50.0, 300.0}) {
      CAPTURE(c);
      CHECK(computeCExp(c) == doctest::Approx(std::exp(c)));
    }
  }

  SUBCASE("Zero where it is not") {
    // zero rather than infinity: arsinhexp's guard never reads the value there, and a masked lane
    // multiplying it by a slip rate of zero would raise the invalid flag on infinity where it
    // stays quiet on zero
    CHECK(computeCExp(1000.0) == 0.0);
    CHECK(computeCExp(logMaxExp<double>()) == 0.0);
    CHECK(computeCExp(200.0F) == 0.0F);
    CHECK(computeCExp(logMaxExp<float>()) == 0.0F);
  }

  SUBCASE("Never infinite, at either precision") {
    for (const double c : {-2000.0, -1000.0, 0.0, 705.0, 1000.0, 2000.0}) {
      CAPTURE(c);
      CHECK(std::isfinite(computeCExp(c)));
      CHECK(std::isfinite(computeCExp(static_cast<float>(c))));
    }
  }

  SUBCASE("A dual argument carries its derivative through") {
    // in a law that folds the state variable, c depends on the slip rate and this exponential is
    // inside the differentiated chain
    const auto result = computeCExp(Dual<double>(3.0, 2.0));
    CHECK(result.value == doctest::Approx(std::exp(3.0)));
    CHECK(result.derivative == doctest::Approx(2.0 * std::exp(3.0)));
    const auto clamped = computeCExp(Dual<double>(1000.0, 2.0));
    CHECK(clamped.value == 0.0);
    CHECK(clamped.derivative == 0.0);
  }
}

// ---------------------------------------------------------------------------
// arsinhexp
// ---------------------------------------------------------------------------

TEST_CASE_TEMPLATE("DR RateAndState arsinhexp over the reachable domain",
                   T,
                   float,
                   double) { // NOLINT
  using namespace rstest;
  // What the friction laws can present: a slip rate between almostZero() and a hundred metres per
  // second, and a state variable over the friction parameter anywhere from deeply negative to a
  // few hundred. Every point of it has to come back with the full precision of T, since mu is
  // this function times the friction parameter and the inversion resolves mu relative to it.
  const double tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-14;
  double worstError = 0.0;
  double worstX = 0.0;
  double worstC = 0.0;
  std::size_t compared = 0;

  for (const double c : sampleCExpLog<T>()) {
    for (const double x : reachableX<T>()) {
      const auto xt = static_cast<T>(x);
      const auto ct = static_cast<T>(c);
      const double target = reference(static_cast<double>(xt), static_cast<double>(ct));
      // only where a T can hold the answer at full precision: below its smallest normal the
      // correctly rounded result is zero or a subnormal, and neither says anything about the
      // formula
      if (!std::isnormal(target) ||
          std::abs(target) < static_cast<double>(std::numeric_limits<T>::min()) ||
          std::abs(target) > static_cast<double>(std::numeric_limits<T>::max())) {
        continue;
      }
      ++compared;
      const double error = relativeError(static_cast<double>(evaluate(xt, ct)), target);
      if (error > worstError) {
        worstError = error;
        worstX = static_cast<double>(xt);
        worstC = static_cast<double>(ct);
      }
    }
  }

  CAPTURE(worstX);
  CAPTURE(worstC);
  CAPTURE(worstError);
  CHECK(compared > 1000);
  CHECK(worstError < tolerance);
}

TEST_CASE("DR RateAndState arsinhexp symmetry and endpoints" *
          doctest::test_suite("dynamicrupture")) {
  using namespace rstest;

  SUBCASE("Odd in x") {
    // the asymptotic branch takes |x| and reapplies the sign, and frexp reads the same exponent
    // for either sign, so that branch is exactly odd; the other inherits whatever symmetry the
    // library asinh has
    for (const double c : sampleCExpLog<double>()) {
      for (const double x : {1e-20, 1e-6, 1.0, 1e6}) {
        CAPTURE(c);
        CAPTURE(x);
        CHECK(evaluate(-x, c) == -evaluate(x, c));
      }
    }
  }

  SUBCASE("Exactly zero at x = 0") {
    // the asymptotic form has a logarithmic singularity at zero that the function it stands in for
    // does not, so this is a guard and not an identity
    for (const double c : sampleCExpLog<double>()) {
      CAPTURE(c);
      CHECK(evaluate(0.0, c) == 0.0);
      CHECK(evaluate(0.0F, static_cast<float>(c)) == 0.0F);
    }
  }

  SUBCASE("Beyond where any exponential is representable") {
    // the only remaining check is the asymptote itself, and with it that nothing overflows on the
    // way: c reaches a few hundred over a friction parameter of 0.01 already
    for (const double c : {800.0, 2000.0, 1e5}) {
      for (const double x : {1e-8, 1.0, 1e8}) {
        CAPTURE(c);
        CAPTURE(x);
        CHECK(evaluate(x, c) == doctest::Approx(c + std::log(2.0 * x)).epsilon(1e-14));
        CHECK(std::isfinite(evaluate(static_cast<float>(x), static_cast<float>(c))));
      }
    }
  }
}

TEST_CASE_TEMPLATE("DR RateAndState arsinhexp is continuous across its switch",
                   T,
                   float,
                   double) { // NOLINT
  using namespace rstest;
  // The guard compares c + max(exponent of x, 0) log 2 against logMaxExp. A wrong branch, or a
  // formula that is not yet asymptotic where it is used, appears as a step at the boundary rather
  // than as a drift, so the boundary is crossed densely in both of the variables that move it.
  const double tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-14;
  const auto limit = static_cast<double>(rs::logMaxExp<T>());

  const auto sweep = [&](double x, double from, double to, double step) {
    double worst = 0.0;
    double worstC = from;
    for (double c = from; c <= to; c += step) {
      const auto xt = static_cast<T>(x);
      const auto ct = static_cast<T>(c);
      const double target = reference(static_cast<double>(xt), static_cast<double>(ct));
      if (!std::isnormal(target)) {
        continue;
      }
      const double error = relativeError(static_cast<double>(evaluate(xt, ct)), target);
      if (error > worst) {
        worst = error;
        worstC = c;
      }
    }
    CAPTURE(x);
    CAPTURE(worstC);
    CAPTURE(worst);
    CHECK(worst < tolerance);
  };

  SUBCASE("In c, at a moderate x") { sweep(1.0, limit - 5.0, limit + 5.0, 1e-3); }
  SUBCASE("In c, at a large x") { sweep(1e6, limit - 20.0, limit + 5.0, 1e-3); }

  SUBCASE("In the exponent of x, at a fixed c") {
    // with c just below the limit, the switch happens as soon as x passes the first binade that
    // lifts the product out of range
    double worst = 0.0;
    double worstX = 0.0;
    for (int i = -400; i <= 400; ++i) {
      const auto xt = static_cast<T>(std::pow(2.0, i / 100.0));
      const auto ct = static_cast<T>(limit - 1.0);
      const double target = reference(static_cast<double>(xt), static_cast<double>(ct));
      if (!std::isnormal(target)) {
        continue;
      }
      const double error = relativeError(static_cast<double>(evaluate(xt, ct)), target);
      if (error > worst) {
        worst = error;
        worstX = static_cast<double>(xt);
      }
    }
    CAPTURE(worstX);
    CAPTURE(worst);
    CHECK(worst < tolerance);
  }
}

TEST_CASE("DR RateAndState arsinhexp argument order is (x, cExpLog, cExp)" *
          doctest::test_suite("dynamicrupture")) {
  using namespace rstest;
  // All three parameters share one type, so a transposed call compiles in silence. The numbers are
  // the ones RateAndStateInitializer forms for a fault at rest.
  constexpr double RsA = 0.008;
  const double x = 5e-11;                    // initialSlipRate * 0.5 / rsSr0
  const double cExpLog = 101.15085092994046; // (f0 + b log(sr0 Psi / sl0)) / a
  const double cExp = rs::computeCExp(cExpLog);

  const double correct = rs::arsinhexp(x, cExpLog, cExp);
  CHECK(correct == doctest::Approx(reference(x, cExpLog)).epsilon(1e-14));
  CHECK(RsA * correct == doctest::Approx(0.625).epsilon(1e-9));

  // the two exponential arguments transposed give a friction coefficient of 1e44 instead of 0.625
  CHECK(relativeError(rs::arsinhexp(x, cExp, cExpLog), correct) > 1e3);
}

TEST_CASE_TEMPLATE("DR RateAndState arsinhexp outside the reachable domain",
                   T,
                   float,
                   double) { // NOLINT
  using namespace rstest;
  // Two bands the header names as inaccurate, both requiring an x the friction laws cannot form.
  // They are pinned rather than left to be discovered: whoever widens the domain will see these
  // fail, and the accuracy they would need costs an exp(c/2) and a second multiplication in the
  // inversion's inner loop.
  const auto limit = static_cast<double>(rs::logMaxExp<T>());
  const auto smallest = static_cast<double>(std::numeric_limits<T>::min());

  SUBCASE("An unrepresentable exp(c) against a subnormal x") {
    // the product is of order one, so the asymptotic branch is entered far from its asymptote
    const auto x = static_cast<T>(smallest / 8.0);
    const auto c = static_cast<T>(limit + 3.0);
    if (x > static_cast<T>(0)) {
      const double target = reference(static_cast<double>(x), static_cast<double>(c));
      CHECK(std::abs(target) > 1e-4);
      CHECK(relativeError(static_cast<double>(evaluate(x, c)), target) > 1e-2);
    }
  }

  SUBCASE("An underflowing exp(c) against a large x") {
    // the mirror image: the precomputed factor is zero and the plain branch returns zero, while
    // the product is still a normal number. exp(c) has to fall below the smallest subnormal for
    // this, which is a decade and a half further out than where it stops being normal.
    const auto c =
        static_cast<T>(std::log(static_cast<double>(std::numeric_limits<T>::denorm_min())) - 5.0);
    // formed through logarithms, since exp(c) underflows a double here as well
    const auto x =
        static_cast<T>(std::exp(std::log(smallest) - static_cast<double>(c) + std::log(100.0)));
    const double target = reference(static_cast<double>(x), static_cast<double>(c));
    REQUIRE(std::isfinite(x));
    REQUIRE(rs::computeCExp(c) == static_cast<T>(0));
    if (std::isnormal(target)) {
      CHECK(evaluate(x, c) == static_cast<T>(0));
      CHECK(relativeError(0.0, target) > 0.5);
    }
  }
}

// ---------------------------------------------------------------------------
// the derivative, which now travels with the value
// ---------------------------------------------------------------------------

TEST_CASE_TEMPLATE("DR RateAndState arsinhexp differentiates by the slip rate",
                   T,
                   float,
                   double) { // NOLINT
  using namespace rstest;
  // The invariant the inversion rests on: once the derivative stops belonging to the value beside
  // it, the Newton step degrades to a secant with no other symptom than a slower iteration.
  const double tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-13;

  SUBCASE("Against the closed form where the product is representable") {
    // exp(c) itself has to be a normal number of T, or the precomputed factor carries the error of
    // a subnormal into everything downstream -- at c = -100 a float has two bits of it left
    const auto limit = static_cast<double>(rs::logMaxExp<T>());
    for (const double c : {-0.9 * limit, -20.0, 0.0, 20.0, 0.7 * limit}) {
      for (const double x : {1e-10, 1e-5, 0.25, 1.0, 7.0, 1e4}) {
        CAPTURE(c);
        CAPTURE(x);
        const auto xt = static_cast<T>(x);
        const auto ct = static_cast<T>(c);
        const auto seeded = Dual<T>(xt, static_cast<T>(1));
        const auto result = rs::arsinhexp(seeded, Dual<T>(ct), rs::computeCExp(Dual<T>(ct)));
        // d/dx asinh(x e^c) = e^c / sqrt(1 + (x e^c)^2)
        const double product = static_cast<double>(xt) * std::exp(static_cast<double>(ct));
        const double target = std::exp(static_cast<double>(ct)) / std::hypot(1.0, product);
        CHECK(static_cast<double>(result.value) ==
              doctest::Approx(reference(static_cast<double>(xt), static_cast<double>(ct))));
        CHECK(relativeError(static_cast<double>(result.derivative), target) < tolerance);
      }
    }
  }

  SUBCASE("Reaching 1 / x where the product does not fit") {
    // past the switch the value is c + log(2x) and the derivative is what that differentiates to,
    // which is also the limit of the closed form -- and the one place where a derivative built on
    // 1 + v*v would have returned zero and asked for an infinite step
    for (const double c : {200.0, 1000.0, 1e4}) {
      for (const double x : {1e-6, 1.0, 1e6}) {
        CAPTURE(c);
        CAPTURE(x);
        const auto xt = static_cast<T>(x);
        const auto ct = static_cast<T>(c);
        const auto result = rs::arsinhexp(
            Dual<T>(xt, static_cast<T>(1)), Dual<T>(ct), rs::computeCExp(Dual<T>(ct)));
        CHECK(relativeError(static_cast<double>(result.derivative), 1.0 / static_cast<double>(xt)) <
              tolerance);
      }
    }
  }

  SUBCASE("A state variable that follows the slip rate contributes its own term") {
    // what a folded law evaluates: c = c0 + k log x, so
    //   d/dx asinh(x e^c) = (1 + k) e^c / sqrt(1 + (x e^c)^2) * ... = (1 + k) / x  for a large
    // product, and the general form below otherwise. A derivative that reached only through x
    // would miss the factor entirely.
    const double slope = 0.35;
    for (const double x : {1e-8, 1.0, 1e4}) {
      CAPTURE(x);
      const auto xt = static_cast<T>(x);
      const auto seeded = Dual<T>(xt, static_cast<T>(1));
      const auto cExpLog =
          Dual<T>(static_cast<T>(20)) + Dual<T>(static_cast<T>(slope)) * log(seeded);
      const auto result = rs::arsinhexp(seeded, cExpLog, rs::computeCExp(cExpLog));

      const double xd = static_cast<double>(xt);
      const double c = 20.0 + slope * std::log(xd);
      const double product = xd * std::exp(c);
      // d/dx [x e^{c(x)}] = e^c (1 + x c'(x)) with c'(x) = slope / x
      const double target = std::exp(c) * (1.0 + slope) / std::hypot(1.0, product);
      CHECK(relativeError(static_cast<double>(result.derivative), target) < tolerance);
    }
  }
}

// ---------------------------------------------------------------------------
// logsinh
// ---------------------------------------------------------------------------

TEST_CASE("DR RateAndState logsinh" * doctest::test_suite("dynamicrupture")) {
  using namespace rstest;

  SUBCASE("Matches log(x sinh(c))") {
    // the production regime: both call sites pass c = |tau / (a p)|, of order one to a hundred, and
    // x = 2 sr0 / V, which is large
    for (const double c : {0.1, 0.5, 1.0, 5.0, 20.0, 50.0, 100.0, 200.0}) {
      for (const double x : {1e-8, 1e-2, 0.5, 2.0, 1e6, 1e12}) {
        CAPTURE(c);
        CAPTURE(x);
        if (std::isfinite(std::sinh(c))) {
          CHECK(relativeError(rs::logsinh(x, c), std::log(x * std::sinh(c))) < 1e-13);
        }
      }
    }
  }

  SUBCASE("Large c asymptotics") {
    // log(x sinh(c)) -> c + log(x/2), where sinh alone would have overflowed past c = 710
    for (const double c : {50.0, 200.0, 700.0, 1500.0}) {
      for (const double x : {1e-8, 2.0, 1e6}) {
        CAPTURE(c);
        CAPTURE(x);
        CHECK(rs::logsinh(x, c) == doctest::Approx(c + std::log(x / 2.0)).epsilon(1e-12));
      }
    }
  }

  SUBCASE("Loses relative accuracy where x sinh(c) approaches one") {
    // |c| + log(...) with the two terms cancelling. At c = 1e-6 and x = 1e6 the result is 1.7e-13
    // and four digits of it are gone. No call site comes near, and the limit should be visible
    // rather than assumed away.
    const double value = rs::logsinh(1e6, 1e-6);
    const double target = std::log(1e6 * std::sinh(1e-6));
    CHECK(relativeError(value, target) > 1e-12);
    CHECK(relativeError(value, target) < 1e-3);
  }

  SUBCASE("Defined exactly where x sinh(c) is positive") {
    CHECK(std::isnan(rs::logsinh(1.0, -1.0)));
    CHECK(relativeError(rs::logsinh(-1.0, -1.0), std::log(-1.0 * std::sinh(-1.0))) < 1e-13);

    // sinh(0) = 0, hence minus infinity rather than a NaN
    const double atZero = rs::logsinh(1.0, 0.0);
    CHECK(std::isinf(atZero));
    CHECK(atZero < 0.0);
  }

  SUBCASE("A dual argument carries its derivative through") {
    // d/dx log(x sinh(c)) = 1/x
    const auto result = rs::logsinh(Dual<double>(3.0, 1.0), Dual<double>(2.0));
    CHECK(result.value == doctest::Approx(std::log(3.0 * std::sinh(2.0))));
    CHECK(result.derivative == doctest::Approx(1.0 / 3.0));
    // and d/dc log(x sinh(c)) = coth(c)
    const auto byC = rs::logsinh(Dual<double>(3.0), Dual<double>(2.0, 1.0));
    CHECK(byC.derivative == doctest::Approx(std::cosh(2.0) / std::sinh(2.0)));
  }
}

// ---------------------------------------------------------------------------
// effectiveNormalStress
// ---------------------------------------------------------------------------

TEST_CASE("DR RateAndState effectiveNormalStress closes the pressurization" *
          doctest::test_suite("dynamicrupture")) {
  using namespace rstest;
  // The claim the closed form rests on is that
  //     sigma = stick / (1 - mu V slope)
  // is the fixed point of
  //     sigma <- stick + slope mu V sigma,
  // which is what an outer iteration would have converged to over several steps. Reaching it in
  // one is the whole point, so the two are compared directly: the iteration is run until it stops
  // moving and has to land on the closed form to the last digit.
  const auto fixedPoint = [](double stick, double slope, double slipRate, double mu) {
    double sigma = stick;
    for (int iteration = 0; iteration < 10000; ++iteration) {
      const double next = stick + slope * mu * slipRate * sigma;
      if (next == sigma) {
        return sigma;
      }
      sigma = next;
    }
    return std::numeric_limits<double>::quiet_NaN();
  };

  SUBCASE("It is the fixed point the outer iteration would have reached") {
    // the slope is negative and the stick compressive; mu V slope stays well inside the unit disc
    // wherever the fault has strength left, which is what makes the iteration converge at all
    for (const double stick : {-1e5, -1e7, -5e7}) {
      for (const double slope : {-1e-12, -1e-9, -1e-8}) {
        for (const double slipRate : {1e-6, 1e-2, 1.0, 10.0}) {
          for (const double mu : {0.1, 0.6, 1.2}) {
            const double contraction = std::abs(slope * mu * slipRate);
            if (contraction > 0.5) {
              continue;
            }
            CAPTURE(stick);
            CAPTURE(slope);
            CAPTURE(slipRate);
            CAPTURE(mu);
            const double closed = rs::effectiveNormalStress(stick, slope, slipRate, mu);
            CHECK(closed == doctest::Approx(fixedPoint(stick, slope, slipRate, mu)).epsilon(1e-14));
          }
        }
      }
    }
  }

  SUBCASE("Heating unloads the fault") {
    // a negative slope lifts the pore pressure, so the effective normal stress has to come out
    // smaller in magnitude than the stick it started from, and monotonically so in the slip rate
    const double stick = -2e7;
    const double slope = -1e-9;
    double previous = std::abs(stick);
    for (const double slipRate : {0.0, 1e-3, 1e-2, 1e-1, 1.0}) {
      CAPTURE(slipRate);
      const double closed = rs::effectiveNormalStress(stick, slope, slipRate, 0.6);
      CHECK(std::abs(closed) <= previous);
      previous = std::abs(closed);
    }
    CHECK(rs::effectiveNormalStress(stick, slope, 0.0, 0.6) == stick);
    CHECK(rs::effectiveNormalStress(stick, 0.0, 1.0, 0.6) == stick);
  }

  SUBCASE("A runaway has no solution on this branch and is clamped") {
    // The divisor is 1 - mu V slope, so with the negative slope of a heating fault it exceeds one
    // for every slip rate and the clamp is a guard rather than a path: no compressive stick can
    // reach it. What does reach it is the opposite sign, where the heating adds to the loading
    // instead of relieving it and the divisor passes through zero.
    const double stick = -2e7;
    CHECK(rs::effectiveNormalStress(stick, 2.0, 1.0, 0.6) == 0.0); // divisor -0.2
    CHECK(rs::effectiveNormalStress(stick, 1.7, 1.0, 0.6) == 0.0); // divisor -0.02
    // just short of it the fault still has a solution, however weak
    const double weak = rs::effectiveNormalStress(stick, 1.6, 1.0, 0.6);
    CHECK(weak < 0.0);
    CHECK(std::abs(weak) > 10.0 * std::abs(stick));
    // a stick that has already gone tensile is clamped whatever the slope does
    CHECK(rs::effectiveNormalStress(1e5, -1e-9, 1.0, 0.6) == 0.0);
    CHECK(rs::effectiveNormalStress(1e5, 2.0, 1.0, 0.6) == 0.0);
    // and with the production sign the result is compressive over the whole range of slip rates
    for (const double slipRate : {0.0, 1e-3, 1.0, 100.0}) {
      CAPTURE(slipRate);
      CHECK(rs::effectiveNormalStress(stick, -1e-9, slipRate, 0.6) < 0.0);
    }
  }

  SUBCASE("The derivative follows the whole closure") {
    // sigma depends on the slip rate directly and through mu, and both reach it through the same
    // divisor. d sigma / dV = sigma^2 slope (mu + V mu') / stick, which a hand-written derivative
    // would have had to carry separately from the value.
    const double stick = -2e7;
    const double slope = -1e-9;
    const double muSlope = 0.02;
    for (const double slipRate : {1e-3, 1.0, 5.0}) {
      CAPTURE(slipRate);
      const Dual<double> trial(slipRate, 1.0);
      const Dual<double> mu = Dual<double>(0.6) + Dual<double>(muSlope) * trial;
      const auto sigma = rs::effectiveNormalStress(Dual<double>(stick), slope, trial, mu);

      const double muValue = 0.6 + muSlope * slipRate;
      const double divisor = 1.0 - muValue * slipRate * slope;
      const double value = stick / divisor;
      // d(divisor)/dV = -slope (mu + V mu'), and d sigma / dV = -stick / divisor^2 * d(divisor)/dV
      const double target = stick / (divisor * divisor) * slope * (muValue + slipRate * muSlope);
      CHECK(sigma.value == doctest::Approx(value));
      CHECK(sigma.derivative == doctest::Approx(target));
    }
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATECOMMON_T_H_
