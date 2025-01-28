// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "tests/TestHelper.h"
#include <array>
#include <iomanip>
#include <iostream>

#include <Eigen/Dense>

#include "Geometry/Refinement/TriangleRefiner.h"

namespace seissol::unit_test {

inline void assertTriangle(seissol::refinement::Triangle& a,
                           std::array<Eigen::Vector2d, 3>& b,
                           double area,
                           double epsilon) {
  REQUIRE(a.area == AbsApprox(area).epsilon(epsilon));
  for (int i = 0; i < 3; i++) {
    REQUIRE(a.x[i][0] == AbsApprox(b[i][0]).epsilon(epsilon));
    REQUIRE(a.x[i][1] == AbsApprox(b[i][1]).epsilon(epsilon));
  }
}

TEST_CASE("Triangle Refiner") {
  constexpr double Epsilon = std::numeric_limits<double>::epsilon();
  seissol::refinement::TriangleRefiner tr;

  SUBCASE("Divide by 1") {
    std::array<std::array<Eigen::Vector2d, 3>, 4> expectedTriangles = {
        {{Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(0.5, 0.0), Eigen::Vector2d(0.0, 0.5)},
         {Eigen::Vector2d(0.5, 0.0), Eigen::Vector2d(1.0, 0.0), Eigen::Vector2d(0.5, 0.5)},
         {Eigen::Vector2d(0.0, 0.5), Eigen::Vector2d(0.5, 0.5), Eigen::Vector2d(0.0, 1.0)},
         {Eigen::Vector2d(0.5, 0.5), Eigen::Vector2d(0.0, 0.5), Eigen::Vector2d(0.5, 0.0)}}};
    const double area = 0.25;

    tr.refine(1);

    for (unsigned i = 0; i < 4; i++) {
      assertTriangle(tr.subTris.at(i), expectedTriangles[i], area, Epsilon);
    }
  }

  SUBCASE("Refine by 3") {
    std::array<std::array<Eigen::Vector2d, 3>, 64> expectedTriangles = {
        {{Eigen::Vector2d(0, 0), Eigen::Vector2d(0.125, 0), Eigen::Vector2d(0, 0.125)},
         {Eigen::Vector2d(0.125, 0), Eigen::Vector2d(0.25, 0), Eigen::Vector2d(0.125, 0.125)},
         {Eigen::Vector2d(0, 0.125), Eigen::Vector2d(0.125, 0.125), Eigen::Vector2d(0, 0.25)},
         {Eigen::Vector2d(0.125, 0.125), Eigen::Vector2d(0, 0.125), Eigen::Vector2d(0.125, 0)},
         {Eigen::Vector2d(0.25, 0), Eigen::Vector2d(0.375, 0), Eigen::Vector2d(0.25, 0.125)},
         {Eigen::Vector2d(0.375, 0), Eigen::Vector2d(0.5, 0), Eigen::Vector2d(0.375, 0.125)},
         {Eigen::Vector2d(0.25, 0.125), Eigen::Vector2d(0.375, 0.125), Eigen::Vector2d(0.25, 0.25)},
         {Eigen::Vector2d(0.375, 0.125), Eigen::Vector2d(0.25, 0.125), Eigen::Vector2d(0.375, 0)},
         {Eigen::Vector2d(0, 0.25), Eigen::Vector2d(0.125, 0.25), Eigen::Vector2d(0, 0.375)},
         {Eigen::Vector2d(0.125, 0.25), Eigen::Vector2d(0.25, 0.25), Eigen::Vector2d(0.125, 0.375)},
         {Eigen::Vector2d(0, 0.375), Eigen::Vector2d(0.125, 0.375), Eigen::Vector2d(0, 0.5)},
         {Eigen::Vector2d(0.125, 0.375), Eigen::Vector2d(0, 0.375), Eigen::Vector2d(0.125, 0.25)},
         {Eigen::Vector2d(0.25, 0.25), Eigen::Vector2d(0.125, 0.25), Eigen::Vector2d(0.25, 0.125)},
         {Eigen::Vector2d(0.125, 0.25), Eigen::Vector2d(0, 0.25), Eigen::Vector2d(0.125, 0.125)},
         {Eigen::Vector2d(0.25, 0.125), Eigen::Vector2d(0.125, 0.125), Eigen::Vector2d(0.25, 0)},
         {Eigen::Vector2d(0.125, 0.125),
          Eigen::Vector2d(0.25, 0.125),
          Eigen::Vector2d(0.125, 0.25)},
         {Eigen::Vector2d(0.5, 0), Eigen::Vector2d(0.625, 0), Eigen::Vector2d(0.5, 0.125)},
         {Eigen::Vector2d(0.625, 0), Eigen::Vector2d(0.75, 0), Eigen::Vector2d(0.625, 0.125)},
         {Eigen::Vector2d(0.5, 0.125), Eigen::Vector2d(0.625, 0.125), Eigen::Vector2d(0.5, 0.25)},
         {Eigen::Vector2d(0.625, 0.125), Eigen::Vector2d(0.5, 0.125), Eigen::Vector2d(0.625, 0)},
         {Eigen::Vector2d(0.75, 0), Eigen::Vector2d(0.875, 0), Eigen::Vector2d(0.75, 0.125)},
         {Eigen::Vector2d(0.875, 0), Eigen::Vector2d(1, 0), Eigen::Vector2d(0.875, 0.125)},
         {Eigen::Vector2d(0.75, 0.125), Eigen::Vector2d(0.875, 0.125), Eigen::Vector2d(0.75, 0.25)},
         {Eigen::Vector2d(0.875, 0.125), Eigen::Vector2d(0.75, 0.125), Eigen::Vector2d(0.875, 0)},
         {Eigen::Vector2d(0.5, 0.25), Eigen::Vector2d(0.625, 0.25), Eigen::Vector2d(0.5, 0.375)},
         {Eigen::Vector2d(0.625, 0.25), Eigen::Vector2d(0.75, 0.25), Eigen::Vector2d(0.625, 0.375)},
         {Eigen::Vector2d(0.5, 0.375), Eigen::Vector2d(0.625, 0.375), Eigen::Vector2d(0.5, 0.5)},
         {Eigen::Vector2d(0.625, 0.375), Eigen::Vector2d(0.5, 0.375), Eigen::Vector2d(0.625, 0.25)},
         {Eigen::Vector2d(0.75, 0.25), Eigen::Vector2d(0.625, 0.25), Eigen::Vector2d(0.75, 0.125)},
         {Eigen::Vector2d(0.625, 0.25), Eigen::Vector2d(0.5, 0.25), Eigen::Vector2d(0.625, 0.125)},
         {Eigen::Vector2d(0.75, 0.125), Eigen::Vector2d(0.625, 0.125), Eigen::Vector2d(0.75, 0)},
         {Eigen::Vector2d(0.625, 0.125),
          Eigen::Vector2d(0.75, 0.125),
          Eigen::Vector2d(0.625, 0.25)},
         {Eigen::Vector2d(0, 0.5), Eigen::Vector2d(0.125, 0.5), Eigen::Vector2d(0, 0.625)},
         {Eigen::Vector2d(0.125, 0.5), Eigen::Vector2d(0.25, 0.5), Eigen::Vector2d(0.125, 0.625)},
         {Eigen::Vector2d(0, 0.625), Eigen::Vector2d(0.125, 0.625), Eigen::Vector2d(0, 0.75)},
         {Eigen::Vector2d(0.125, 0.625), Eigen::Vector2d(0, 0.625), Eigen::Vector2d(0.125, 0.5)},
         {Eigen::Vector2d(0.25, 0.5), Eigen::Vector2d(0.375, 0.5), Eigen::Vector2d(0.25, 0.625)},
         {Eigen::Vector2d(0.375, 0.5), Eigen::Vector2d(0.5, 0.5), Eigen::Vector2d(0.375, 0.625)},
         {Eigen::Vector2d(0.25, 0.625), Eigen::Vector2d(0.375, 0.625), Eigen::Vector2d(0.25, 0.75)},
         {Eigen::Vector2d(0.375, 0.625), Eigen::Vector2d(0.25, 0.625), Eigen::Vector2d(0.375, 0.5)},
         {Eigen::Vector2d(0, 0.75), Eigen::Vector2d(0.125, 0.75), Eigen::Vector2d(0, 0.875)},
         {Eigen::Vector2d(0.125, 0.75), Eigen::Vector2d(0.25, 0.75), Eigen::Vector2d(0.125, 0.875)},
         {Eigen::Vector2d(0, 0.875), Eigen::Vector2d(0.125, 0.875), Eigen::Vector2d(0, 1)},
         {Eigen::Vector2d(0.125, 0.875), Eigen::Vector2d(0, 0.875), Eigen::Vector2d(0.125, 0.75)},
         {Eigen::Vector2d(0.25, 0.75), Eigen::Vector2d(0.125, 0.75), Eigen::Vector2d(0.25, 0.625)},
         {Eigen::Vector2d(0.125, 0.75), Eigen::Vector2d(0, 0.75), Eigen::Vector2d(0.125, 0.625)},
         {Eigen::Vector2d(0.25, 0.625), Eigen::Vector2d(0.125, 0.625), Eigen::Vector2d(0.25, 0.5)},
         {Eigen::Vector2d(0.125, 0.625),
          Eigen::Vector2d(0.25, 0.625),
          Eigen::Vector2d(0.125, 0.75)},
         {Eigen::Vector2d(0.5, 0.5), Eigen::Vector2d(0.375, 0.5), Eigen::Vector2d(0.5, 0.375)},
         {Eigen::Vector2d(0.375, 0.5), Eigen::Vector2d(0.25, 0.5), Eigen::Vector2d(0.375, 0.375)},
         {Eigen::Vector2d(0.5, 0.375), Eigen::Vector2d(0.375, 0.375), Eigen::Vector2d(0.5, 0.25)},
         {Eigen::Vector2d(0.375, 0.375), Eigen::Vector2d(0.5, 0.375), Eigen::Vector2d(0.375, 0.5)},
         {Eigen::Vector2d(0.25, 0.5), Eigen::Vector2d(0.125, 0.5), Eigen::Vector2d(0.25, 0.375)},
         {Eigen::Vector2d(0.125, 0.5), Eigen::Vector2d(0, 0.5), Eigen::Vector2d(0.125, 0.375)},
         {Eigen::Vector2d(0.25, 0.375), Eigen::Vector2d(0.125, 0.375), Eigen::Vector2d(0.25, 0.25)},
         {Eigen::Vector2d(0.125, 0.375), Eigen::Vector2d(0.25, 0.375), Eigen::Vector2d(0.125, 0.5)},
         {Eigen::Vector2d(0.5, 0.25), Eigen::Vector2d(0.375, 0.25), Eigen::Vector2d(0.5, 0.125)},
         {Eigen::Vector2d(0.375, 0.25), Eigen::Vector2d(0.25, 0.25), Eigen::Vector2d(0.375, 0.125)},
         {Eigen::Vector2d(0.5, 0.125), Eigen::Vector2d(0.375, 0.125), Eigen::Vector2d(0.5, 0)},
         {Eigen::Vector2d(0.375, 0.125), Eigen::Vector2d(0.5, 0.125), Eigen::Vector2d(0.375, 0.25)},
         {Eigen::Vector2d(0.25, 0.25), Eigen::Vector2d(0.375, 0.25), Eigen::Vector2d(0.25, 0.375)},
         {Eigen::Vector2d(0.375, 0.25), Eigen::Vector2d(0.5, 0.25), Eigen::Vector2d(0.375, 0.375)},
         {Eigen::Vector2d(0.25, 0.375), Eigen::Vector2d(0.375, 0.375), Eigen::Vector2d(0.25, 0.5)},
         {Eigen::Vector2d(0.375, 0.375),
          Eigen::Vector2d(0.25, 0.375),
          Eigen::Vector2d(0.375, 0.25)}}};
    const double area = 0.015625;

    tr.refine(3);

    for (unsigned i = 0; i < 64; i++) {
      assertTriangle(tr.subTris.at(i), expectedTriangles[i], area, Epsilon);
    }
  }
}

} // namespace seissol::unit_test
