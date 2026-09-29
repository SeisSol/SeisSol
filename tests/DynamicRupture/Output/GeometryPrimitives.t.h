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

TEST_CASE("ExtTriangle construction and access" * doctest::test_suite("dynamicrupture")) {
  const CoordinateT p0{0.0, 0.0, 0.0};
  const CoordinateT p1{1.0, 0.0, 0.0};
  const CoordinateT p2{0.0, 1.0, 0.0};
  ExtTriangle tri(p0, p1, p2);

  CHECK(tri.point(0)[0] == doctest::Approx(0.0));
  CHECK(tri.point(1)[0] == doctest::Approx(1.0));
  CHECK(tri.point(2)[1] == doctest::Approx(1.0));
  CHECK(ExtTriangle::size() == 3);
}

} // namespace seissol::unit_test
