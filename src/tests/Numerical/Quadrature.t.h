// SPDX-FileCopyrightText: 2017-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Numerical/Quadrature.h"

namespace seissol::unit_test {

TEST_CASE("Test generation of Gauss Jacobi quadrature points") {
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

TEST_CASE("test triangle quadrature") {
  constexpr auto Epsilon = std::numeric_limits<double>::epsilon() * 10;

  double points[4][2];
  double weights[4];
  seissol::quadrature::TriangleQuadrature(points, weights, 2);
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

} // namespace seissol::unit_test
