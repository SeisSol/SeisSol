// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_DUAL_T_H_
#define SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_DUAL_T_H_

#include <doctest.h>

#include "DynamicRupture/FrictionLaws/Dual.h"
#include "Kernels/Precision.h"
#include "TestHelper.h"

#include <cmath>
#include <limits>
#include <type_traits>
#include <vector>

namespace seissol::unit_test {

namespace dualtest {

using seissol::dr::friction_law::derivativeOf;
using seissol::dr::friction_law::Dual;
using seissol::dr::friction_law::dualCast;
using seissol::dr::friction_law::valueOf;

/// What a derivative assembled from a value and a rule may lose against a reference evaluated in
/// double: the elementary function itself, one or two products, and one division.
template <typename T>
constexpr double tolerance() {
  return std::is_same_v<T, float> ? 2e-6 : 1e-13;
}

/**
 * A purely relative comparison. doctest's own Approx carries an absolute floor of one epsilon
 * scale, which would accept any value at all once the reference drops below the tolerance -- and
 * several derivatives here are 1e-16 or smaller, which is exactly where they have to be checked.
 */
inline AbsApprox relative(double reference, double tol) {
  return AbsApprox(reference)
      .epsilon(reference == 0.0 ? std::numeric_limits<double>::min() : 0.0)
      .delta(tol);
}

/**
 * Seeds the differentiation variable at each point of a list and compares both halves of the
 * result against closed forms evaluated in double.
 *
 * The reference is evaluated at the point the function receives, not at the point the list names.
 * Near a singularity the two are not interchangeable: the derivative of asin at 0.999999 changes
 * by seven parts in a thousand over one single-precision step, so a reference taken at the
 * unrounded argument would measure how well 0.999999 fits in a float and nothing else.
 */
template <typename T, typename Function, typename Value, typename Derivative>
void unary(const std::vector<double>& points,
           Function function,
           Value value,
           Derivative derivative) {
  for (const double x : points) {
    CAPTURE(x);
    const auto seed = static_cast<T>(x);
    const auto received = static_cast<double>(seed);
    const auto result = function(Dual<T>(seed, static_cast<T>(1)));
    CHECK(static_cast<double>(result.value) == relative(value(received), tolerance<T>()));
    CHECK(static_cast<double>(result.derivative) == relative(derivative(received), tolerance<T>()));
  }
}

/// Differentiates a two-argument function in one slot at a time, holding the other constant.
template <typename T, typename Function, typename First, typename Second>
void binary(const std::vector<std::pair<double, double>>& points,
            Function function,
            First first,
            Second second) {
  for (const auto& [x, y] : points) {
    CAPTURE(x);
    CAPTURE(y);
    const auto firstSeed = static_cast<T>(x);
    const auto secondSeed = static_cast<T>(y);
    const auto a = static_cast<double>(firstSeed);
    const auto b = static_cast<double>(secondSeed);
    const auto byFirst = function(Dual<T>(firstSeed, static_cast<T>(1)), Dual<T>(secondSeed));
    const auto bySecond = function(Dual<T>(firstSeed), Dual<T>(secondSeed, static_cast<T>(1)));
    CHECK(static_cast<double>(byFirst.derivative) == relative(first(a, b), tolerance<T>()));
    CHECK(static_cast<double>(bySecond.derivative) == relative(second(a, b), tolerance<T>()));
  }
}

} // namespace dualtest

// ---------------------------------------------------------------------------
// plumbing
// ---------------------------------------------------------------------------

TEST_CASE("DR Dual seeding" * doctest::test_suite("dynamicrupture")) {
  using namespace dualtest;

  SUBCASE("A constant carries no dependence on the variable") {
    CHECK(Dual<double>(7.0).value == 7.0);
    CHECK(Dual<double>(7.0).derivative == 0.0);
  }

  SUBCASE("valueOf and derivativeOf read either kind of scalar") {
    CHECK(valueOf(3.5) == 3.5);
    CHECK(derivativeOf(3.5) == 0.0);
    CHECK(valueOf(Dual<double>(3.5, 1.5)) == 3.5);
    CHECK(derivativeOf(Dual<double>(3.5, 1.5)) == 1.5);
  }

  SUBCASE("dualCast changes the precision and keeps the derivative") {
    const auto narrowed = dualCast<float>(Dual<double>(1.5, 2.5));
    CHECK(narrowed.value == 1.5F);
    CHECK(narrowed.derivative == 2.5F);
    CHECK(dualCast<float>(1.5) == 1.5F);
  }
}

TEST_CASE("DR Dual arithmetic" * doctest::test_suite("dynamicrupture")) {
  using namespace dualtest;
  const Dual<double> a(2.0, 3.0);
  const Dual<double> b(5.0, 7.0);

  SUBCASE("The four operations") {
    CHECK((a + b).value == 7.0);
    CHECK((a + b).derivative == 10.0);
    CHECK((a - b).value == -3.0);
    CHECK((a - b).derivative == -4.0);
    CHECK((a * b).value == 10.0);
    CHECK((a * b).derivative == 3.0 * 5.0 + 2.0 * 7.0);
    CHECK((a / b).value == doctest::Approx(0.4));
    CHECK((a / b).derivative == doctest::Approx((3.0 * 5.0 - 2.0 * 7.0) / 25.0));
    CHECK((+a).derivative == 3.0);
    CHECK((-a).derivative == -3.0);
  }

  SUBCASE("A plain scalar on either side") {
    // deduction rejects a mixed pair without these, which is what the friction laws would hit
    // first were they to drop a cast
    CHECK((a + 5.0).value == 7.0);
    CHECK((a + 5.0).derivative == 3.0);
    CHECK((5.0 + a).value == 7.0);
    CHECK((5.0 + a).derivative == 3.0);
    CHECK((a - 5.0).value == -3.0);
    CHECK((a - 5.0).derivative == 3.0);
    CHECK((5.0 - a).value == 3.0);
    CHECK((5.0 - a).derivative == -3.0);
    CHECK((a * 5.0).derivative == 15.0);
    CHECK((5.0 * a).derivative == 15.0);
    CHECK((a / 4.0).derivative == 0.75);
    CHECK((4.0 / a).value == 2.0);
    CHECK((4.0 / a).derivative == -3.0);
  }

  SUBCASE("Compound assignment matches the binary operator") {
    Dual<double> sum(2.0, 3.0);
    sum += b;
    CHECK(sum.value == (a + b).value);
    CHECK(sum.derivative == (a + b).derivative);

    Dual<double> difference(2.0, 3.0);
    difference -= b;
    CHECK(difference.derivative == (a - b).derivative);

    Dual<double> product(2.0, 3.0);
    product *= b;
    CHECK(product.derivative == (a * b).derivative);

    Dual<double> quotient(2.0, 3.0);
    quotient /= b;
    CHECK(quotient.value == (a / b).value);
    CHECK(quotient.derivative == (a / b).derivative);
  }
}

TEST_CASE("DR Dual classification" * doctest::test_suite("dynamicrupture")) {
  using namespace dualtest;
  using seissol::dr::friction_law::isfinite;
  using seissol::dr::friction_law::isinf;
  using seissol::dr::friction_law::isnan;
  const double infinity = std::numeric_limits<double>::infinity();

  // a finite value beside an infinite derivative is what a step that has left the domain looks
  // like, and the assertions in the friction laws are there to catch it
  CHECK(isfinite(Dual<double>(1.0, 2.0)));
  CHECK_FALSE(isfinite(Dual<double>(1.0, infinity)));
  CHECK_FALSE(isfinite(Dual<double>(infinity, 2.0)));
  CHECK(isnan(Dual<double>(1.0, std::nan(""))));
  CHECK_FALSE(isnan(Dual<double>(1.0, 2.0)));
  CHECK(isinf(Dual<double>(infinity, 0.0)));
  CHECK(isinf(Dual<double>(0.0, -infinity)));
}

// ---------------------------------------------------------------------------
// elementary functions, against closed forms
// ---------------------------------------------------------------------------

TEST_CASE_TEMPLATE("DR Dual exponentials and logarithms", T, float, double) { // NOLINT
  using namespace dualtest;
  constexpr double Ln2 = 0.69314718055994530942;
  constexpr double Ln10 = 2.30258509299404568402;

  unary<T>(
      {-3.0, -0.5, 0.0, 0.5, 3.0, 20.0},
      [](auto d) { return exp(d); },
      [](double x) { return std::exp(x); },
      [](double x) { return std::exp(x); });
  unary<T>(
      {-3.0, -0.5, 0.0, 0.5, 3.0, 20.0},
      [](auto d) { return exp2(d); },
      [](double x) { return std::exp2(x); },
      [Ln2](double x) { return std::exp2(x) * Ln2; });
  unary<T>(
      {-3.0, -1e-6, 0.0, 1e-6, 3.0},
      [](auto d) { return expm1(d); },
      [](double x) { return std::expm1(x); },
      [](double x) { return std::exp(x); });
  unary<T>(
      {1e-6, 0.5, 1.0, 2.0, 1e6},
      [](auto d) { return log(d); },
      [](double x) { return std::log(x); },
      [](double x) { return 1.0 / x; });
  unary<T>(
      {1e-6, 0.5, 1.0, 2.0, 1e6},
      [](auto d) { return log2(d); },
      [](double x) { return std::log2(x); },
      [Ln2](double x) { return 1.0 / (x * Ln2); });
  unary<T>(
      {1e-6, 0.5, 1.0, 2.0, 1e6},
      [](auto d) { return log10(d); },
      [](double x) { return std::log10(x); },
      [Ln10](double x) { return 1.0 / (x * Ln10); });
  unary<T>(
      {-0.9, -1e-6, 0.0, 1e-6, 3.0, 1e6},
      [](auto d) { return log1p(d); },
      [](double x) { return std::log1p(x); },
      [](double x) { return 1.0 / (1.0 + x); });
}

TEST_CASE_TEMPLATE("DR Dual powers and roots", T, float, double) { // NOLINT
  using namespace dualtest;

  unary<T>(
      {1e-6, 0.25, 1.0, 4.0, 1e6},
      [](auto d) { return sqrt(d); },
      [](double x) { return std::sqrt(x); },
      [](double x) { return 0.5 / std::sqrt(x); });
  unary<T>(
      {-8.0, -1e-6, 1e-6, 1.0, 8.0, 1e6},
      [](auto d) { return cbrt(d); },
      [](double x) { return std::cbrt(x); },
      [](double x) { return 1.0 / (3.0 * std::cbrt(x) * std::cbrt(x)); });

  // a constant exponent, which is the shape the friction laws use
  for (const double exponent : {-1.5, 0.125, 3.0}) {
    for (const double x : {0.5, 2.0, 1e3}) {
      CAPTURE(exponent);
      CAPTURE(x);
      const auto result =
          pow(Dual<T>(static_cast<T>(x), static_cast<T>(1)), static_cast<T>(exponent));
      CHECK(static_cast<double>(result.value) == relative(std::pow(x, exponent), tolerance<T>()));
      CHECK(static_cast<double>(result.derivative) ==
            relative(exponent * std::pow(x, exponent - 1.0), tolerance<T>()));
    }
  }

  // a constant base
  for (const double base : {0.5, 2.0, 10.0}) {
    for (const double exponent : {-1.5, 0.125, 3.0}) {
      CAPTURE(base);
      CAPTURE(exponent);
      const auto result =
          pow(static_cast<T>(base), Dual<T>(static_cast<T>(exponent), static_cast<T>(1)));
      CHECK(static_cast<double>(result.derivative) ==
            relative(std::pow(base, exponent) * std::log(base), tolerance<T>()));
    }
  }

  // both, which needs the logarithm of the base as well
  binary<T>(
      {{2.0, 3.0}, {2.0, -1.5}, {0.5, 0.25}, {1e3, 0.3}},
      [](auto x, auto y) { return pow(x, y); },
      [](double a, double b) { return b * std::pow(a, b - 1.0); },
      [](double a, double b) { return std::pow(a, b) * std::log(a); });

  binary<T>(
      {{3.0, 4.0}, {-3.0, 4.0}, {3.0, -4.0}, {0.5, 0.25}},
      [](auto x, auto y) { return hypot(x, y); },
      [](double a, double b) { return a / std::hypot(a, b); },
      [](double a, double b) { return b / std::hypot(a, b); });
}

TEST_CASE_TEMPLATE("DR Dual trigonometric functions", T, float, double) { // NOLINT
  using namespace dualtest;

  unary<T>(
      {-3.0, -0.5, 0.0, 0.5, 1.5, 3.0},
      [](auto d) { return sin(d); },
      [](double x) { return std::sin(x); },
      [](double x) { return std::cos(x); });
  unary<T>(
      {-3.0, -0.5, 0.0, 0.5, 1.5, 3.0},
      [](auto d) { return cos(d); },
      [](double x) { return std::cos(x); },
      [](double x) { return -std::sin(x); });
  unary<T>(
      {-1.0, -0.5, 0.0, 0.5, 1.0, 1.5},
      [](auto d) { return tan(d); },
      [](double x) { return std::tan(x); },
      [](double x) { return 1.0 / (std::cos(x) * std::cos(x)); });
  unary<T>(
      {-0.999, -0.5, 0.0, 0.5, 0.999},
      [](auto d) { return asin(d); },
      [](double x) { return std::asin(x); },
      [](double x) { return 1.0 / std::sqrt((1.0 - x) * (1.0 + x)); });
  unary<T>(
      {-0.999, -0.5, 0.0, 0.5, 0.999},
      [](auto d) { return acos(d); },
      [](double x) { return std::acos(x); },
      [](double x) { return -1.0 / std::sqrt((1.0 - x) * (1.0 + x)); });
  unary<T>(
      {-1e3, -1.0, 0.0, 1.0, 1e3},
      [](auto d) { return atan(d); },
      [](double x) { return std::atan(x); },
      [](double x) { return 1.0 / (1.0 + x * x); });

  // atan2 picks the quadrant from both signs, and both derivatives change sign with it
  binary<T>(
      {{3.0, 4.0}, {-3.0, 4.0}, {3.0, -4.0}, {-3.0, -4.0}, {0.5, 0.25}},
      [](auto y, auto x) { return atan2(y, x); },
      [](double y, double x) { return x / (x * x + y * y); },
      [](double y, double x) { return -y / (x * x + y * y); });
}

TEST_CASE_TEMPLATE("DR Dual hyperbolic functions", T, float, double) { // NOLINT
  using namespace dualtest;

  unary<T>(
      {-5.0, -0.5, 0.0, 0.5, 5.0},
      [](auto d) { return sinh(d); },
      [](double x) { return std::sinh(x); },
      [](double x) { return std::cosh(x); });
  unary<T>(
      {-5.0, -0.5, 0.0, 0.5, 5.0},
      [](auto d) { return cosh(d); },
      [](double x) { return std::cosh(x); },
      [](double x) { return std::sinh(x); });
  unary<T>(
      {-5.0, -0.5, 0.0, 0.5, 5.0},
      [](auto d) { return tanh(d); },
      [](double x) { return std::tanh(x); },
      [](double x) { return 1.0 / (std::cosh(x) * std::cosh(x)); });
  unary<T>(
      {-1e3, -1.0, 0.0, 1.0, 1e3},
      [](auto d) { return asinh(d); },
      [](double x) { return std::asinh(x); },
      [](double x) { return 1.0 / std::hypot(1.0, x); });
  unary<T>(
      {1.0001, 1.5, 2.0, 1e3},
      [](auto d) { return acosh(d); },
      [](double x) { return std::acosh(x); },
      [](double x) { return 1.0 / (std::sqrt(x - 1.0) * std::sqrt(x + 1.0)); });
  unary<T>(
      {-0.999, -0.5, 0.0, 0.5, 0.999},
      [](auto d) { return atanh(d); },
      [](double x) { return std::atanh(x); },
      [](double x) { return 1.0 / ((1.0 - x) * (1.0 + x)); });
}

TEST_CASE_TEMPLATE("DR Dual error functions", T, float, double) { // NOLINT
  using namespace dualtest;
  constexpr double TwoOverSqrtPi = 1.12837916709551257390;

  unary<T>(
      {-3.0, -0.5, 0.0, 0.5, 3.0},
      [](auto d) { return erf(d); },
      [](double x) { return std::erf(x); },
      [TwoOverSqrtPi](double x) { return TwoOverSqrtPi * std::exp(-x * x); });
  unary<T>(
      {-3.0, -0.5, 0.0, 0.5, 3.0},
      [](auto d) { return erfc(d); },
      [](double x) { return std::erfc(x); },
      [TwoOverSqrtPi](double x) { return -TwoOverSqrtPi * std::exp(-x * x); });
}

// ---------------------------------------------------------------------------
// the piecewise functions, where the derivative is a choice
// ---------------------------------------------------------------------------

TEST_CASE("DR Dual selection and sign" * doctest::test_suite("dynamicrupture")) {
  using namespace dualtest;
  using seissol::dr::friction_law::abs;
  using seissol::dr::friction_law::copysign;
  using seissol::dr::friction_law::fabs;
  using seissol::dr::friction_law::fdim;
  using seissol::dr::friction_law::fmax;
  using seissol::dr::friction_law::fmin;

  SUBCASE("The selected operand keeps its own derivative") {
    // which is what a friction law clamping its steady-state friction coefficient needs: below
    // the clamp the coefficient is frozen and contributes nothing to the slope
    CHECK(fmax(Dual<double>(2.0, 5.0), Dual<double>(3.0, 7.0)).value == 3.0);
    CHECK(fmax(Dual<double>(2.0, 5.0), Dual<double>(3.0, 7.0)).derivative == 7.0);
    CHECK(fmin(Dual<double>(2.0, 5.0), Dual<double>(3.0, 7.0)).value == 2.0);
    CHECK(fmin(Dual<double>(2.0, 5.0), Dual<double>(3.0, 7.0)).derivative == 5.0);
  }

  SUBCASE("A tie resolves to the first operand") {
    CHECK(fmax(Dual<double>(2.0, 5.0), Dual<double>(2.0, 7.0)).derivative == 5.0);
    CHECK(fmin(Dual<double>(2.0, 5.0), Dual<double>(2.0, 7.0)).derivative == 5.0);
  }

  SUBCASE("abs is the derivative from the right at the kink") {
    CHECK(abs(Dual<double>(-3.0, 2.0)).value == 3.0);
    CHECK(abs(Dual<double>(-3.0, 2.0)).derivative == -2.0);
    CHECK(abs(Dual<double>(3.0, 2.0)).derivative == 2.0);
    CHECK(abs(Dual<double>(0.0, 1.0)).derivative == 1.0);
    CHECK(fabs(Dual<double>(-3.0, 2.0)).derivative == -2.0);
  }

  SUBCASE("copysign carries the product of the two signs") {
    CHECK(copysign(Dual<double>(3.0, 2.0), Dual<double>(-1.0)).value == -3.0);
    CHECK(copysign(Dual<double>(3.0, 2.0), Dual<double>(-1.0)).derivative == -2.0);
    CHECK(copysign(Dual<double>(-3.0, 2.0), Dual<double>(-1.0)).derivative == 2.0);
    CHECK(copysign(Dual<double>(-3.0, 2.0), Dual<double>(1.0)).derivative == -2.0);
  }

  SUBCASE("fdim closes, with its derivative") {
    CHECK(fdim(Dual<double>(5.0, 3.0), Dual<double>(2.0, 1.0)).value == 3.0);
    CHECK(fdim(Dual<double>(5.0, 3.0), Dual<double>(2.0, 1.0)).derivative == 2.0);
    CHECK(fdim(Dual<double>(1.0, 3.0), Dual<double>(2.0, 1.0)).value == 0.0);
    CHECK(fdim(Dual<double>(1.0, 3.0), Dual<double>(2.0, 1.0)).derivative == 0.0);
  }
}

TEST_CASE("DR Dual piecewise constant functions" * doctest::test_suite("dynamicrupture")) {
  using namespace dualtest;
  using seissol::dr::friction_law::ceil;
  using seissol::dr::friction_law::floor;
  using seissol::dr::friction_law::round;
  using seissol::dr::friction_law::trunc;

  // constant between their steps, and the derivative says so
  CHECK(floor(Dual<double>(-2.5, 1.0)).value == -3.0);
  CHECK(floor(Dual<double>(-2.5, 1.0)).derivative == 0.0);
  CHECK(ceil(Dual<double>(-2.5, 1.0)).value == -2.0);
  CHECK(ceil(Dual<double>(-2.5, 1.0)).derivative == 0.0);
  CHECK(round(Dual<double>(2.4, 1.0)).value == 2.0);
  CHECK(round(Dual<double>(2.4, 1.0)).derivative == 0.0);
  CHECK(trunc(Dual<double>(-2.5, 1.0)).value == -2.0);
  CHECK(trunc(Dual<double>(-2.5, 1.0)).derivative == 0.0);
}

TEST_CASE("DR Dual composites" * doctest::test_suite("dynamicrupture")) {
  using namespace dualtest;
  using seissol::dr::friction_law::fma;
  using seissol::dr::friction_law::ldexp;

  CHECK(fma(Dual<double>(2.0, 1.0), Dual<double>(3.0, 1.0), Dual<double>(4.0, 1.0)).value == 10.0);
  CHECK(fma(Dual<double>(2.0, 1.0), Dual<double>(3.0, 1.0), Dual<double>(4.0, 1.0)).derivative ==
        1.0 * 3.0 + 2.0 * 1.0 + 1.0);

  // scaling by a power of two is exact, so the derivative scales with it and nothing else changes
  CHECK(ldexp(Dual<double>(3.0, 2.0), 4).value == 48.0);
  CHECK(ldexp(Dual<double>(3.0, 2.0), 4).derivative == 32.0);
  CHECK(ldexp(Dual<double>(3.0, 2.0), -2).derivative == 0.5);
}

// ---------------------------------------------------------------------------
// the properties that matter to the inversion
// ---------------------------------------------------------------------------

TEST_CASE("DR Dual the chain rule composes" * doctest::test_suite("dynamicrupture")) {
  using namespace dualtest;
  // an identity built from two functions has to come back with derivative one exactly, which no
  // single rule can fake
  for (const double x : {0.5, 2.0, 30.0}) {
    CAPTURE(x);
    CHECK(log(exp(Dual<double>(x, 1.0))).derivative == doctest::Approx(1.0));
    CHECK(exp(log(Dual<double>(x, 1.0))).derivative == doctest::Approx(1.0));
    CHECK(sqrt(Dual<double>(x, 1.0) * Dual<double>(x, 1.0)).derivative == doctest::Approx(1.0));
    CHECK(sinh(asinh(Dual<double>(x, 1.0))).derivative == doctest::Approx(1.0));
  }

  SUBCASE("A product of three factors gets all three terms") {
    const Dual<double> x(2.0, 1.0);
    const auto cube = x * x * x;
    CHECK(cube.value == 8.0);
    CHECK(cube.derivative == doctest::Approx(3.0 * 4.0));
  }
}

TEST_CASE_TEMPLATE("DR Dual derivatives keep the range of their value",
                   T,
                   float,
                   double) { // NOLINT
  using namespace dualtest;
  // Each of these is a point where the textbook expression for the derivative loses the answer
  // while the function itself is untroubled. A zero where the derivative should be finite is the
  // worst failure available to a Newton iteration: it asks for an infinite step.
  const auto edge = static_cast<double>(std::is_same_v<T, float> ? 1e20F : 1e150);

  SUBCASE("asinh at the edge of the range") {
    // 1 + v*v overflows from v = 1.8e19 in single and 1.3e154 in double precision; asinh itself
    // does not overflow at any representable argument, and arsinhexp evaluates it out here
    unary<T>(
        {edge, -edge},
        [](auto d) { return asinh(d); },
        [](double x) { return std::asinh(x); },
        [](double x) { return 1.0 / std::hypot(1.0, x); });
  }

  SUBCASE("acosh at a large argument") {
    // v*v - 1 overflows at the same place
    unary<T>(
        {edge},
        [](auto d) { return acosh(d); },
        [](double x) { return std::acosh(x); },
        [](double x) { return 1.0 / (std::sqrt(x - 1.0) * std::sqrt(x + 1.0)); });
  }

  SUBCASE("tanh once the tangent has reached one") {
    // 1 - tanh(v)^2 is zero in single precision from v = 9 and in double from v = 19, while
    // 1 / cosh(v)^2 stays exact until cosh itself squares out of range
    unary<T>(
        {12.0, -12.0},
        [](auto d) { return tanh(d); },
        [](double x) { return std::tanh(x); },
        [](double x) { return 1.0 / (std::cosh(x) * std::cosh(x)); });
  }

  SUBCASE("hypot and atan2 with a leg at the top of the range") {
    // the numerator a*a' would overflow before the quotient does
    binary<T>(
        {{edge, 2.0 * edge}, {edge, 0.0}},
        [](auto x, auto y) { return hypot(x, y); },
        [](double a, double b) { return a / std::hypot(a, b); },
        [](double a, double b) { return b / std::hypot(a, b); });
    binary<T>(
        {{edge, 2.0 * edge}},
        [](auto y, auto x) { return atan2(y, x); },
        [](double y, double x) { return (x / std::hypot(y, x)) / std::hypot(y, x); },
        [](double y, double x) { return -(y / std::hypot(y, x)) / std::hypot(y, x); });
  }

  SUBCASE("asin and acos at the ends of their interval") {
    // (1-v)(1+v) keeps the digits that 1 - v*v drops, and the derivative diverges here, so the
    // digits are the whole answer
    const double near = 1.0 - (std::is_same_v<T, float> ? 1e-6 : 1e-12);
    unary<T>(
        {near, -near},
        [](auto d) { return asin(d); },
        [](double x) { return std::asin(x); },
        [](double x) { return 1.0 / std::sqrt((1.0 - x) * (1.0 + x)); });
  }
}

TEST_CASE("DR Dual in working precision" * doctest::test_suite("dynamicrupture")) {
  using namespace dualtest;
  // the friction laws instantiate the type at `real`, whichever that is in this build
  unary<real>(
      {0.25, 1.0, 3.0},
      [](auto d) { return exp(d); },
      [](double x) { return std::exp(x); },
      [](double x) { return std::exp(x); });
  unary<real>(
      {0.25, 1.0, 3.0},
      [](auto d) { return asinh(d); },
      [](double x) { return std::asinh(x); },
      [](double x) { return 1.0 / std::hypot(1.0, x); });
  CHECK(valueOf(Dual<real>(static_cast<real>(2), static_cast<real>(1))) == static_cast<real>(2));
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_DUAL_T_H_
