// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "doctest.h"

#include "Config.h"
#include "GeneratedCode/init.h"
#include "IO/Instance/Geometry/Points.h"
#include "Kernels/Precision.h"
#include "TestHelper.h"

#include <limits>

namespace seissol::unit_test {

TEST_CASE("IO/Points") {
  const auto pointsCompare = [](auto pointsView, const auto& generated) {
    REQUIRE(pointsView.shape(1) == generated.size());
    REQUIRE(pointsView.shape(0) == generated[0].size());
    for (std::size_t pointId = 0; pointId < generated.size(); ++pointId) {
      for (std::size_t coord = 0; coord < generated[pointId].size(); ++coord) {
        const auto ref = pointsView.isInRange(coord, pointId) ? pointsView(coord, pointId) : 0;
        REQUIRE(ref == AbsApprox(generated[pointId][coord])
                           .epsilon(std::numeric_limits<real>::epsilon())
                           .delta(std::numeric_limits<real>::epsilon()));
      }
    }
  };

  SUBCASE("Triangle (2D)") {
    using Vtk2d = init::vtk2d<Config>;
    using seissol::io::instance::geometry::pointsTriangle;
    pointsCompare(Vtk2d::view<1>::create(Vtk2d::Values1), pointsTriangle(1));
    pointsCompare(Vtk2d::view<2>::create(Vtk2d::Values2), pointsTriangle(2));
    pointsCompare(Vtk2d::view<3>::create(Vtk2d::Values3), pointsTriangle(3));
    pointsCompare(Vtk2d::view<4>::create(Vtk2d::Values4), pointsTriangle(4));
    pointsCompare(Vtk2d::view<5>::create(Vtk2d::Values5), pointsTriangle(5));
    pointsCompare(Vtk2d::view<6>::create(Vtk2d::Values6), pointsTriangle(6));
    pointsCompare(Vtk2d::view<7>::create(Vtk2d::Values7), pointsTriangle(7));
    pointsCompare(Vtk2d::view<8>::create(Vtk2d::Values8), pointsTriangle(8));
  }

  SUBCASE("Tetrahedron (3D)") {
    using Vtk3d = init::vtk3d<Config>;
    using seissol::io::instance::geometry::pointsTetrahedron;
    pointsCompare(Vtk3d::view<1>::create(Vtk3d::Values1), pointsTetrahedron(1));
    pointsCompare(Vtk3d::view<2>::create(Vtk3d::Values2), pointsTetrahedron(2));
    pointsCompare(Vtk3d::view<3>::create(Vtk3d::Values3), pointsTetrahedron(3));
    pointsCompare(Vtk3d::view<4>::create(Vtk3d::Values4), pointsTetrahedron(4));
    pointsCompare(Vtk3d::view<5>::create(Vtk3d::Values5), pointsTetrahedron(5));
    pointsCompare(Vtk3d::view<6>::create(Vtk3d::Values6), pointsTetrahedron(6));
    pointsCompare(Vtk3d::view<7>::create(Vtk3d::Values7), pointsTetrahedron(7));
    pointsCompare(Vtk3d::view<8>::create(Vtk3d::Values8), pointsTetrahedron(8));
  }
}

} // namespace seissol::unit_test
