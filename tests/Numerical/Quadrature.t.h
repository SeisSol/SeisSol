// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Numerical/Quadrature.h"
#include "Numerical/SymmetricTetrahedronRule.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("Test generation of Gauss Jacobi quadrature points" * doctest::test_suite("numerical")) {
  constexpr auto Epsilon = std::numeric_limits<double>::epsilon() * 10;
  double points[5];
  double weights[5];
  seissol::quadrature::GaussJacobi(points, weights, 5, 1, 3);
  // Compare to Maple reference solution
  CHECK(points[0] == AbsApprox(0.86698568210542769702).epsilon(Epsilon));
  CHECK(points[1] == AbsApprox(0.57652877512667440772).epsilon(Epsilon));
  CHECK(points[2] == AbsApprox(0.17976783188823737401).epsilon(Epsilon));
  CHECK(points[3] == AbsApprox(-.25499675973326581341).epsilon(Epsilon));
  CHECK(points[4] == AbsApprox(-.65399981510135937963).epsilon(Epsilon));
  CHECK(weights[0] == AbsApprox(0.18915446768616357329).epsilon(Epsilon));
  CHECK(weights[1] == AbsApprox(0.58714974961811369751).epsilon(Epsilon));
  CHECK(weights[2] == AbsApprox(0.57657004957734461768).epsilon(Epsilon));
  CHECK(weights[3] == AbsApprox(0.22255926867518051648).epsilon(Epsilon));
  CHECK(weights[4] == AbsApprox(0.024566464443197594119).epsilon(Epsilon));
}

TEST_CASE("test triangle quadrature" * doctest::test_suite("numerical")) {
  constexpr auto Epsilon = std::numeric_limits<double>::epsilon() * 10;

  std::array<std::array<double, 2>, 4> points{};
  std::array<double, 4> weights{};
  seissol::quadrature::TriangleQuadrature(points.data(), weights.data(), 2);
  // Compare to Maple reference solution
  CHECK(points[0][0] == AbsApprox(0.64494897427831780982).epsilon(Epsilon));
  CHECK(points[1][0] == AbsApprox(0.64494897427831780982).epsilon(Epsilon));
  CHECK(points[2][0] == AbsApprox(0.15505102572168219018).epsilon(Epsilon));
  CHECK(points[3][0] == AbsApprox(0.15505102572168219018).epsilon(Epsilon));
  CHECK(points[0][1] == AbsApprox(0.28001991549907407200).epsilon(Epsilon));
  CHECK(points[1][1] == AbsApprox(0.075031110222608118175).epsilon(Epsilon));
  CHECK(points[2][1] == AbsApprox(0.66639024601470138669).epsilon(Epsilon));
  CHECK(points[3][1] == AbsApprox(0.17855872826361642311).epsilon(Epsilon));
  CHECK(weights[0] == AbsApprox(0.090979309128011415315).epsilon(Epsilon));
  CHECK(weights[1] == AbsApprox(0.090979309128011415315).epsilon(Epsilon));
  CHECK(weights[2] == AbsApprox(0.15902069087198858472).epsilon(Epsilon));
  CHECK(weights[3] == AbsApprox(0.15902069087198858472).epsilon(Epsilon));
}

TEST_CASE("Symmetric tetrahedron quadrature rules" * doctest::test_suite("numerical")) {
  const auto factorial = [](std::size_t n) {
    double value = 1;
    for (std::size_t i = 2; i <= n; ++i) {
      value *= static_cast<double>(i);
    }
    return value;
  };
  const auto barycentric = [](const std::array<double, 3>& x) {
    return std::array<double, 4>{1 - x[0] - x[1] - x[2], x[0], x[1], x[2]};
  };

  for (std::size_t degree = 0; degree <= seissol::quadrature::MaxSymmetricTetrahedronDegree;
       ++degree) {
    const auto [points, weights] = seissol::quadrature::symmetricTetrahedronRule(degree);
    REQUIRE(points.size() == seissol::quadrature::symmetricTetrahedronRuleSize(degree));
    REQUIRE(weights.size() == points.size());

    for (std::size_t i = 0; i < points.size(); ++i) {
      const auto lambda = barycentric(points[i]);
      CHECK(weights[i] > 0);
      CHECK(*std::min_element(lambda.begin(), lambda.end()) > 1e-3);
    }

    // exact for x^a y^b z^c with a + b + c <= degree
    for (std::size_t a = 0; a <= degree; ++a) {
      for (std::size_t b = 0; a + b <= degree; ++b) {
        for (std::size_t c = 0; a + b + c <= degree; ++c) {
          const auto exact = factorial(a) * factorial(b) * factorial(c) / factorial(a + b + c + 3);
          double value = 0;
          for (std::size_t i = 0; i < points.size(); ++i) {
            value += weights[i] * std::pow(points[i][0], a) * std::pow(points[i][1], b) *
                     std::pow(points[i][2], c);
          }
          CHECK(std::abs(value - exact) <= 1e-12 * exact);
        }
      }
    }

    // every permutation of the vertices maps the rule onto itself
    std::array<std::size_t, 4> permutation{0, 1, 2, 3};
    do {
      for (std::size_t i = 0; i < points.size(); ++i) {
        const auto lambda = barycentric(points[i]);
        bool found = false;
        for (std::size_t j = 0; j < points.size() && !found; ++j) {
          const auto other = barycentric(points[j]);
          bool same = std::abs(weights[i] - weights[j]) < 1e-14;
          for (std::size_t k = 0; k < 4; ++k) {
            same = same && std::abs(lambda[permutation[k]] - other[k]) < 1e-12;
          }
          found = same;
        }
        CHECK(found);
      }
    } while (std::next_permutation(permutation.begin(), permutation.end()));
  }
  CHECK(seissol::quadrature::symmetricTetrahedronRuleSize(
            seissol::quadrature::MaxSymmetricTetrahedronDegree + 1) == 0);
}

} // namespace seissol::unit_test
