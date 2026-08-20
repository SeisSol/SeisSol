// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "TestHelper.h"

#include <Eigen/Dense>
#include <Geometry/CellTransform.h>
#include <Geometry/MeshDefinition.h>

namespace seissol::unit_test {

TEST_CASE("Affine transform function conformacy") {
  CoordinateT v1{1, 2, 3};
  CoordinateT v2{1, 2, 3};
  CoordinateT v3{1, 2, 3};
  CoordinateT v4{1, 2, 3};

  const auto transform = geometry::AffineTransform({v1, v2, v3, v4});

  // transform.refToSpace(TODO);
}

} // namespace seissol::unit_test
