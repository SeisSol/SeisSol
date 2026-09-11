// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "DynamicRupture/Output/Geometry.h"

namespace seissol::unit_test {
using seissol::dr::ExtTriangle;
using seissol::dr::ExtVrtxCoords;

TEST_CASE("ExtVrtxCoords initializer list" * doctest::test_suite("dynamicrupture")) {
  ExtVrtxCoords p{1.0, 2.0, 3.0};
  CHECK(p[0] == doctest::Approx(1.0));
  CHECK(p[1] == doctest::Approx(2.0));
  CHECK(p[2] == doctest::Approx(3.0));
}

TEST_CASE("ExtVrtxCoords default is zero" * doctest::test_suite("dynamicrupture")) {
  ExtVrtxCoords p;
  CHECK(p[0] == doctest::Approx(0.0));
  CHECK(p[1] == doctest::Approx(0.0));
  CHECK(p[2] == doctest::Approx(0.0));
}

TEST_CASE("ExtVrtxCoords getAsEigen3LibVector" * doctest::test_suite("dynamicrupture")) {
  ExtVrtxCoords p{4.0, 5.0, 6.0};
  auto v = p.getAsEigen3LibVector();
  CHECK(v(0) == doctest::Approx(4.0));
  CHECK(v(1) == doctest::Approx(5.0));
  CHECK(v(2) == doctest::Approx(6.0));
}

TEST_CASE("ExtVrtxCoords size" * doctest::test_suite("dynamicrupture")) {
  CHECK(ExtVrtxCoords::size() == 3);
}

TEST_CASE("ExtTriangle construction and access" * doctest::test_suite("dynamicrupture")) {
  ExtVrtxCoords p0{0.0, 0.0, 0.0};
  ExtVrtxCoords p1{1.0, 0.0, 0.0};
  ExtVrtxCoords p2{0.0, 1.0, 0.0};
  ExtTriangle tri(p0, p1, p2);

  CHECK(tri.point(0)[0] == doctest::Approx(0.0));
  CHECK(tri.point(1)[0] == doctest::Approx(1.0));
  CHECK(tri.point(2)[1] == doctest::Approx(1.0));
  CHECK(ExtTriangle::size() == 3);
}

} // namespace seissol::unit_test
