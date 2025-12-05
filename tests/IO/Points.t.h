// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "doctest.h"

#include "GeneratedCode/init.h"
#include "IO/Instance/Geometry/Points.h"
#include "Kernels/Precision.h"
#include "TestHelper.h"

namespace seissol::unit_test {

TEST_CASE("IO/Points") {
  const auto pointsCompare = [](auto pointsView, const auto& generated) {
    REQUIRE(pointsView.shape(1) == generated.size());
    REQUIRE(pointsView.shape(0) == generated[0].size());
    for (std::size_t pointId = 0; pointId < generated.size(); ++pointId) {
      for (std::size_t coord = 0; coord < generated[pointId].size(); ++coord) {
        const auto ref = pointsView.isInRange(coord, pointId) ? pointsView(coord, pointId) : 0;
        REQUIRE(ref == AbsApprox(generated[pointId][coord]));
      }
    }
  };

  SUBCASE("Triangle (2D)") {
    using init::vtk2d;
    using io::instance::geometry::pointsTriangle;
    pointsCompare(vtk2d::view<1>::create(const_cast<real*>(vtk2d::Values1)), pointsTriangle(1));
    pointsCompare(vtk2d::view<2>::create(const_cast<real*>(vtk2d::Values2)), pointsTriangle(2));
    pointsCompare(vtk2d::view<3>::create(const_cast<real*>(vtk2d::Values3)), pointsTriangle(3));
    pointsCompare(vtk2d::view<4>::create(const_cast<real*>(vtk2d::Values4)), pointsTriangle(4));
    pointsCompare(vtk2d::view<5>::create(const_cast<real*>(vtk2d::Values5)), pointsTriangle(5));
    pointsCompare(vtk2d::view<6>::create(const_cast<real*>(vtk2d::Values6)), pointsTriangle(6));
    pointsCompare(vtk2d::view<7>::create(const_cast<real*>(vtk2d::Values7)), pointsTriangle(7));
    pointsCompare(vtk2d::view<8>::create(const_cast<real*>(vtk2d::Values8)), pointsTriangle(8));
  }

  SUBCASE("Tetrahedron (3D)") {
    using init::vtk3d;
    using io::instance::geometry::pointsTetrahedron;
    pointsCompare(vtk3d::view<1>::create(const_cast<real*>(vtk3d::Values1)), pointsTetrahedron(1));
    pointsCompare(vtk3d::view<2>::create(const_cast<real*>(vtk3d::Values2)), pointsTetrahedron(2));
    pointsCompare(vtk3d::view<3>::create(const_cast<real*>(vtk3d::Values3)), pointsTetrahedron(3));
    pointsCompare(vtk3d::view<4>::create(const_cast<real*>(vtk3d::Values4)), pointsTetrahedron(4));
    pointsCompare(vtk3d::view<5>::create(const_cast<real*>(vtk3d::Values5)), pointsTetrahedron(5));
    pointsCompare(vtk3d::view<6>::create(const_cast<real*>(vtk3d::Values6)), pointsTetrahedron(6));
    pointsCompare(vtk3d::view<7>::create(const_cast<real*>(vtk3d::Values7)), pointsTetrahedron(7));
    pointsCompare(vtk3d::view<8>::create(const_cast<real*>(vtk3d::Values8)), pointsTetrahedron(8));
  }
}

} // namespace seissol::unit_test
