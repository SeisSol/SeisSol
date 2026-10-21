// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Kernels/Precision.h"
#include "Numerical/Functions.h"
#include "Numerical/Projection.h"
#include "Numerical/Quadrature.h"
#include "Numerical/StableSum.h"
#include "TestHelper.h"

namespace seissol::unit_test {

TEST_CASE("Test Jacobi polynomials" * doctest::test_suite("numerical")) {
  constexpr double Epsilon = 10 * std::numeric_limits<double>::epsilon();
  // Compare to Maple reference solution
  CHECK(seissol::functions::JacobiP(0, 1, 0, 0.5) == AbsApprox(1.0).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiP(1, 0, 0, -0.3) ==
        AbsApprox(-.3000000000000000).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiP(1, 13, 4, -1.0) ==
        AbsApprox(-5.000000000000000).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiP(10, 0, 4, -1.0) == AbsApprox(1001.).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiP(10, 0, 4, -0.3) ==
        AbsApprox(-1.870307686685156).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiP(10, 0, 4, 0.0) ==
        AbsApprox(0.8671875000000000).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiP(10, 0, 4, 0.3) ==
        AbsApprox(-.3404032083257812).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiP(10, 0, 4, 1.0) == AbsApprox(1.).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiP(2, 1, 4, -1.0) == AbsApprox(15.).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiP(3, 2, 3, -0.3) ==
        AbsApprox(1.249375000000000).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiP(4, 3, 2, 0.0) ==
        AbsApprox(0.9375000000000002).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiP(5, 4, 1, 0.3) ==
        AbsApprox(-1.982012812500000).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiP(6, 5, 0, 1.0) == AbsApprox(462.0000000000000).epsilon(Epsilon));
}

TEST_CASE("Test first derivative of Jacobi polynomials" * doctest::test_suite("numerical")) {
  constexpr double Epsilon = 100 * std::numeric_limits<real>::epsilon();
  // Compare to Maple reference solution
  CHECK(seissol::functions::JacobiPDerivative(0, 1, 0, 0.5) == AbsApprox(0.0).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiPDerivative(1, 0, 0, -0.3) ==
        AbsApprox(1.000000000000000).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiPDerivative(1, 13, 4, -0.9) ==
        AbsApprox(9.500000000000000).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiPDerivative(7, 0, 4, -0.3) ==
        AbsApprox(-9.122014125000001).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiPDerivative(7, 0, 4, 0.0) ==
        AbsApprox(7.8750000000000).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiPDerivative(7, 0, 4, 0.3) ==
        AbsApprox(-6.067144125000002).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiPDerivative(7, 0, 4, 0.9) ==
        AbsApprox(-2.71497712500000).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiPDerivative(2, 1, 4, -0.9) ==
        AbsApprox(-22.20000000000000).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiPDerivative(3, 2, 3, -0.3) ==
        AbsApprox(3.318750000000000).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiPDerivative(4, 3, 2, 0.0) ==
        AbsApprox(-3.749999999999998).epsilon(Epsilon));
  CHECK(seissol::functions::JacobiPDerivative(5, 4, 1, 0.3) ==
        AbsApprox(-15.37232812500000).epsilon(Epsilon));
}

TEST_CASE("Test Dubiner polynomials" * doctest::test_suite("numerical")) {
  constexpr double Epsilon = 10 * std::numeric_limits<double>::epsilon();
  const std::vector<std::array<double, 7>> tests = {
      {0, 0, 0, 0.25, 0.25, 0.25, 1},
      {1, 0, 0, 0.25, 0.25, 0.25, 0.},
      {0, 1, 0, 0.25, 0.25, 0.25, 0.},
      {0, 0, 1, 0.25, 0.25, 0.25, 0.},
      {2, 0, 0, 0.25, 0.25, 0.25, -0.1250},
      {1, 1, 0, 0.25, 0.25, 0.25, 0.},
      {0, 2, 0, 0.25, 0.25, 0.25, -0.3125},
      {1, 0, 1, 0.25, 0.25, 0.25, 0.},
      {0, 1, 1, 0.25, 0.25, 0.25, 0.},
      {0, 0, 2, 0.25, 0.25, 0.25, -0.5625},
      {3, 0, 0, 0.25, 0.25, 0.25, -0.},
      {2, 1, 0, 0.25, 0.25, 0.25, -0.125000},
      {1, 2, 0, 0.25, 0.25, 0.25, -0.},
      {0, 3, 0, 0.25, 0.25, 0.25, 0.125000},
      {2, 0, 1, 0.25, 0.25, 0.25, -0.125000},
      {1, 1, 1, 0.25, 0.25, 0.25, 0.},
      {0, 2, 1, 0.25, 0.25, 0.25, -0.31250000000000000000},
      {1, 0, 2, 0.25, 0.25, 0.25, -0.},
      {0, 1, 2, 0.25, 0.25, 0.25, -0.},
      {0, 0, 3, 0.25, 0.25, 0.25, 0.437500},
      {0, 0, 0, 0., 0.5, 0.5, 1.},
      {1, 0, 0, 0., 0.5, 0.5, 0.},
      {0, 1, 0, 0., 0.5, 0.5, 1.0},
      {0, 0, 1, 0., 0.5, 0.5, 1.0},
      {2, 0, 0, 0., 0.5, 0.5, 0.},
      {1, 1, 0, 0., 0.5, 0.5, 0.},
      {0, 2, 0, 0., 0.5, 0.5, 0.75},
      {1, 0, 1, 0., 0.5, 0.5, 0.},
      {0, 1, 1, 0., 0.5, 0.5, 2.00},
      {0, 0, 2, 0., 0.5, 0.5, -0.25},
      {3, 0, 0, 0., 0.5, 0.5, 0.},
      {2, 1, 0, 0., 0.5, 0.5, 0.},
      {1, 2, 0, 0., 0.5, 0.5, 0.},
      {0, 3, 0, 0., 0.5, 0.5, 0.500},
      {2, 0, 1, 0., 0.5, 0.5, 0.},
      {1, 1, 1, 0., 0.5, 0.5, 0.},
      {0, 2, 1, 0., 0.5, 0.5, 2.2500000000000000000},
      {1, 0, 2, 0., 0.5, 0.5, 0.},
      {0, 1, 2, 0., 0.5, 0.5, 1.0},
      {0, 0, 3, 0., 0.5, 0.5, -0.750},
      {0, 0, 0, 0., 1.0, 0., 1.},
      {1, 0, 0, 0., 1.0, 0., 0.},
      {0, 1, 0, 0., 1.0, 0., 2.0},
      {0, 0, 1, 0., 1.0, 0., -1.},
      {2, 0, 0, 0., 1.0, 0., 0.},
      {1, 1, 0, 0., 1.0, 0., 0.},
      {0, 2, 0, 0., 1.0, 0., 3.00},
      {1, 0, 1, 0., 1.0, 0., -0.},
      {0, 1, 1, 0., 1.0, 0., -2.0},
      {0, 0, 2, 0., 1.0, 0., 1.},
      {3, 0, 0, 0., 1.0, 0., 0.},
      {2, 1, 0, 0., 1.0, 0., 0.},
      {1, 2, 0, 0., 1.0, 0., 0.},
      {0, 3, 0, 0., 1.0, 0., 4.000},
      {2, 0, 1, 0., 1.0, 0., -0.},
      {1, 1, 1, 0., 1.0, 0., -0.},
      {0, 2, 1, 0., 1.0, 0., -3.0000000000000000000},
      {1, 0, 2, 0., 1.0, 0., 0.},
      {0, 1, 2, 0., 1.0, 0., 2.0},
      {0, 0, 3, 0., 1.0, 0., -1.}};
  for (const auto& t : tests) {
    CHECK(
        seissol::functions::TetraDubinerP(
            {static_cast<unsigned>(t[0]), static_cast<unsigned>(t[1]), static_cast<unsigned>(t[2])},
            {t[3], t[4], t[5]}) == AbsApprox(t[6]).epsilon(Epsilon));
  }
}

TEST_CASE("Bernstein polynomials" * doctest::test_suite("numerical")) {
  constexpr double Epsilon = 10 * std::numeric_limits<double>::epsilon();
  // 4! / (1! 1! 1! 1!) / 4^4 at the barycentre
  CHECK(functions::BernsteinP<3>(4, {1, 1, 1}, {0.25, 0.25, 0.25}) ==
        AbsApprox(24.0 / 256.0).epsilon(Epsilon));
  // 3! / (0! 2! 1!) 0.2^2 0.3, with nothing left for the vertex at the origin
  CHECK(functions::BernsteinP<2>(3, {2, 1}, {0.2, 0.3}) == AbsApprox(0.036).epsilon(Epsilon));
  // 2! / (1! 1!) 0.5^2 on the interval
  CHECK(functions::BernsteinP<1>(2, {1}, {0.5}) == AbsApprox(0.5).epsilon(Epsilon));
  // the polynomial of a vertex is one there and zero at the other vertices
  CHECK(functions::BernsteinP<3>(3, {0, 3, 0}, {0.0, 1.0, 0.0}) == AbsApprox(1.0).epsilon(Epsilon));
  CHECK(functions::BernsteinP<3>(3, {0, 0, 0}, {0.0, 1.0, 0.0}) == AbsApprox(0.0).epsilon(Epsilon));
  CHECK(functions::BernsteinP<3>(3, {0, 0, 0}, {0.0, 0.0, 0.0}) == AbsApprox(1.0).epsilon(Epsilon));

  // a partition of unity of nonnegative functions on the element, so the gradients sum to zero;
  // and the gradients are the ones of the values. The points are binary fractions, so that the
  // barycentric coordinates of the ones on the boundary are zero, not a rounding below it.
  const std::vector<std::array<double, 3>> points = {{0.25, 0.25, 0.25},
                                                     {0.125, 0.25, 0.375},
                                                     {0.5, 0.125, 0.0625},
                                                     {0.0, 0.0, 1.0},
                                                     {0.25, 0.0, 0.75}};
  constexpr double Step = 1e-6;
  for (unsigned degree = 0; degree <= 6; ++degree) {
    const auto indices = numerical::projection::modalIndices<3>(degree + 1);
    for (const auto& point : points) {
      double sum = 0.0;
      std::array<double, 3> gradientSum{};
      for (const auto& index : indices) {
        const double value = functions::BernsteinP<3>(degree, index, point);
        CHECK(value >= 0.0);
        sum += value;
        const auto gradient = functions::gradBernsteinP<3>(degree, index, point);
        for (std::size_t d = 0; d < 3; ++d) {
          gradientSum[d] += gradient[d];
          auto forward = point;
          auto backward = point;
          forward[d] += Step;
          backward[d] -= Step;
          const double difference = (functions::BernsteinP<3>(degree, index, forward) -
                                     functions::BernsteinP<3>(degree, index, backward)) /
                                    (2.0 * Step);
          CHECK(gradient[d] == AbsApprox(difference).epsilon(1e-7));
        }
      }
      CHECK(sum == AbsApprox(1.0).epsilon(100 * Epsilon));
      for (std::size_t d = 0; d < 3; ++d) {
        CHECK(gradientSum[d] == AbsApprox(0.0).epsilon(1000 * Epsilon));
      }
    }
  }
}

TEST_CASE("Bases on the tetrahedron and the triangle" * doctest::test_suite("numerical")) {
  constexpr double Epsilon = 1e-11;
  // the normalized Dubiner functions are orthonormal; the rule of 5 points per dimension
  // integrates a product of two functions of degree 4 exactly
  const auto tetrahedron = seissol::quadrature::simplexRule<3>(5);
  const auto indices3 = numerical::projection::modalIndices<3>(5);
  for (const auto& left : indices3) {
    for (const auto& right : indices3) {
      double product = 0.0;
      for (std::size_t q = 0; q < tetrahedron.first.size(); ++q) {
        product +=
            tetrahedron.second[q] *
            functions::BasisP<BasisType::DubinerNormalized, 3>(4, left, tetrahedron.first[q]) *
            functions::BasisP<BasisType::DubinerNormalized, 3>(4, right, tetrahedron.first[q]);
      }
      CHECK(product == AbsApprox(left == right ? 1.0 : 0.0).epsilon(Epsilon));
    }
  }
  const auto triangle = seissol::quadrature::simplexRule<2>(5);
  const auto indices2 = numerical::projection::modalIndices<2>(5);
  for (const auto& left : indices2) {
    for (const auto& right : indices2) {
      double product = 0.0;
      for (std::size_t q = 0; q < triangle.first.size(); ++q) {
        product += triangle.second[q] *
                   functions::BasisP<BasisType::DubinerNormalized, 2>(4, left, triangle.first[q]) *
                   functions::BasisP<BasisType::DubinerNormalized, 2>(4, right, triangle.first[q]);
      }
      CHECK(product == AbsApprox(left == right ? 1.0 : 0.0).epsilon(Epsilon));
    }
  }
  // the Dubiner basis is DubinerP itself, and a normalized gradient is the Dubiner one divided
  // by the norm of the function
  const std::array<double, 3> point{0.1, 0.2, 0.3};
  for (const auto& index : indices3) {
    CHECK(functions::BasisP<BasisType::Dubiner, 3>(4, index, point) ==
          functions::DubinerP<3>(index, point));
    const auto scaled = functions::gradBasisP<BasisType::DubinerNormalized, 3>(4, index, point);
    const auto plain = functions::gradDubinerP<3>(index, point);
    for (std::size_t d = 0; d < 3; ++d) {
      CHECK(scaled[d] * std::sqrt(functions::DubinerNormSquared<3>(index)) ==
            AbsApprox(plain[d]).epsilon(Epsilon));
    }
  }
}

TEST_CASE("Legendre Polynomials" * doctest::test_suite("numerical")) {
  constexpr double Epsilon = 10 * std::numeric_limits<double>::epsilon();
  constexpr double EpsilonRelaxed = 1e8 * std::numeric_limits<double>::epsilon();

  // compare some common contraints first
  CHECK(seissol::functions::shiftedLegendre(0, 0.5, 1) == AbsApprox(0.0).epsilon(Epsilon));

  CHECK(seissol::functions::shiftedLegendre(0, 0.5, 2) == AbsApprox(0.0).epsilon(Epsilon));
  CHECK(seissol::functions::shiftedLegendre(1, 0.5, 2) == AbsApprox(0.0).epsilon(Epsilon));

  CHECK(seissol::functions::shiftedLegendre(0, 0, 0) ==
        AbsApprox(seissol::functions::shiftedLegendre(0, 1.0, 0)).epsilon(Epsilon));
  CHECK(seissol::functions::shiftedLegendre(1, 0, 1) ==
        AbsApprox(seissol::functions::shiftedLegendre(1, 1.0, 1)).epsilon(Epsilon));
  CHECK(seissol::functions::shiftedLegendre(2, 0, 2) ==
        AbsApprox(seissol::functions::shiftedLegendre(2, 1.0, 2)).epsilon(Epsilon));

  // compare derivatives
  constexpr auto Spacing = 1000;
  for (unsigned n = 0; n < 10; ++n) {
    // test pointwise
    for (int i = 0; i <= Spacing; ++i) {
      const auto x = i / static_cast<double>(Spacing);
      CHECK(seissol::functions::shiftedLegendre(n, x, 0) ==
            AbsApprox(functions::DubinerP<1>({n}, {x})).epsilon(EpsilonRelaxed));
      CHECK(seissol::functions::shiftedLegendre(n, x, 1) ==
            AbsApprox(functions::gradDubinerP<1>({n}, {x})[0]).epsilon(EpsilonRelaxed));
    }

    // test integral
    numerical::StableAccumulator<double> acc;
    for (int i = 0; i < Spacing; ++i) {
      const auto x = i / static_cast<double>(Spacing);
      const auto xp = (i + 1) / static_cast<double>(Spacing);

      // reference via trapezoid
      acc += (seissol::functions::shiftedLegendre(n, x, 0) +
              seissol::functions::shiftedLegendre(n, xp, 0)) /
             static_cast<double>(2 * Spacing);

      // (note: _very_ high error tolerance; due to the trapezoidal sum)
      CHECK(seissol::functions::shiftedLegendre(n, xp, -1) -
                seissol::functions::shiftedLegendre(n, 0, -1) ==
            AbsApprox(acc.result()).epsilon(0.001));
    }
  }
}

} // namespace seissol::unit_test
