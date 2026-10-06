// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Constants.h"
#include "Geometry/CellGeometry.h"
#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"
#include "Geometry/IsoparametricTransform.h"
#include "MockReader.h"
#include "TestHelper.h"

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <vector>

namespace seissol::unit_test::cellgeometry {

using seissol::geometry::AffineFaceTransform;
using seissol::geometry::AffineTransform;
using seissol::geometry::IsoparametricTransform;
using VectorT = seissol::geometry::CellTransform::VectorEigenT;
using FaceVectorT = seissol::geometry::FaceTransform::FaceVectorT;

TEST_CASE("Cell geometry of a mesh" * doctest::test_suite("geometry")) {
  constexpr double Epsilon = 1e-12;
  const std::array<Eigen::Vector3d, Cell::NumVertices> vertices{Eigen::Vector3d(0.1, 0.2, -0.3),
                                                                Eigen::Vector3d(2.1, 0.0, 0.4),
                                                                Eigen::Vector3d(-0.2, 1.7, 0.1),
                                                                Eigen::Vector3d(0.3, 0.1, 2.2)};
  MockReader mesh(vertices);
  const std::vector<VectorT> points{
      VectorT(0.1, 0.2, 0.3), VectorT(0.25, 0.25, 0.25), VectorT(0.7, 0.1, 0.05)};

  SUBCASE("A straight-sided mesh gives the affine transforms") {
    REQUIRE(mesh.geometryOrder() == 1);
    const auto affine = AffineTransform::fromMeshCell(0, mesh);
    const auto transform = seissol::geometry::cellTransformOf(0, mesh);
    for (const auto& point : points) {
      REQUIRE((transform->refToSpace(point) - affine.refToSpace(point)).norm() ==
              AbsApprox(0.0).epsilon(Epsilon));
    }
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const auto face = seissol::geometry::faceTransformOf(0, side, mesh);
      const auto expected = AffineFaceTransform::fromMeshCell(0, side, mesh);
      const auto center = FaceVectorT(Face::ReferenceBarycenter.data());
      REQUIRE((face->normal(center) - expected.normal(center)).norm() ==
              AbsApprox(0.0).epsilon(Epsilon));
    }
  }

  // bend the cell: every edge midpoint moves off the straight edge
  std::array<VectorT, Cell::NumVertices> corners{};
  for (std::size_t i = 0; i < Cell::NumVertices; ++i) {
    corners[i] = vertices[i];
  }
  std::array<VectorT, 6> midpoints{};
  constexpr std::array<std::array<std::size_t, 2>, 6> Edges{
      {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
  for (std::size_t edge = 0; edge < Edges.size(); ++edge) {
    midpoints[edge] = 0.5 * (corners[Edges[edge][0]] + corners[Edges[edge][1]]) +
                      0.05 * VectorT(1.0 + edge, -0.5 * edge, 0.3);
  }
  const auto curved = IsoparametricTransform::fromEdgeMidpoints(corners, midpoints);
  std::vector<CoordinateT> nodes;
  for (const auto& node : curved.nodes()) {
    nodes.push_back({node(0), node(1), node(2)});
  }

  SUBCASE("A curved mesh gives the transform through its nodes") {
    mesh.setCurvedGeometry(2, nodes);
    REQUIRE(mesh.geometryOrder() == 2);
    REQUIRE(mesh.cellNodes(0).size() == nodes.size());
    const auto transform = seissol::geometry::cellTransformOf(0, mesh);
    for (const auto& point : points) {
      REQUIRE((transform->refToSpace(point) - curved.refToSpace(point)).norm() ==
              AbsApprox(0.0).epsilon(Epsilon));
      REQUIRE((transform->refToSpaceJacobian(point) - curved.refToSpaceJacobian(point)).norm() ==
              AbsApprox(0.0).epsilon(Epsilon));
    }
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const auto face = seissol::geometry::faceTransformOf(0, side, mesh);
      const FaceVectorT point(0.2, 0.3);
      REQUIRE((face->refToSpace(point) - curved.refToSpace(face->refToCell(point))).norm() ==
              AbsApprox(0.0).epsilon(Epsilon));
    }
  }

  SUBCASE("Moving the mesh moves the nodes of its curved cells along") {
    mesh.setCurvedGeometry(2, nodes);
    const Eigen::Vector3d displacement(1.0, -2.0, 0.5);
    mesh.displaceMesh(displacement);
    Eigen::Matrix3d scaling;
    scaling << 2.0, 0.1, 0.0, 0.0, 1.5, 0.0, 0.2, 0.0, 0.5;
    mesh.scaleMesh(scaling);
    const auto transform = seissol::geometry::cellTransformOf(0, mesh);
    for (const auto& point : points) {
      const VectorT expected = scaling * (curved.refToSpace(point) + displacement);
      REQUIRE((transform->refToSpace(point) - expected).norm() == AbsApprox(0.0).epsilon(Epsilon));
    }
  }
}

} // namespace seissol::unit_test::cellgeometry
