// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Constants.h"
#include "Geometry/CellTransform.h"
#include "Geometry/FaceTransform.h"
#include "Geometry/MeshTools.h"
#include "TestHelper.h"

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <random>
#include <vector>

namespace seissol::unit_test {

namespace {
using seissol::geometry::AffineFaceTransform;
using seissol::geometry::AffineTransform;
using seissol::geometry::FaceOrientation;
using seissol::geometry::ReferenceFaceMap;

using CellVectorT = seissol::geometry::CellTransform::VectorEigenT;
using FaceVectorT = seissol::geometry::FaceTransform::FaceVectorT;

constexpr std::array<FaceOrientation, 4> AllOrientations{FaceOrientation::Local,
                                                         FaceOrientation::Rotate0,
                                                         FaceOrientation::Rotate1,
                                                         FaceOrientation::Rotate2};

/// the reference tetrahedron, which the reference face maps have to be consistent with
auto referenceCellVertices() -> std::array<CellVectorT, Cell::NumVertices> {
  return {CellVectorT(0, 0, 0), CellVectorT(1, 0, 0), CellVectorT(0, 1, 0), CellVectorT(0, 0, 1)};
}

auto samplePointsOnReferenceFace(std::mt19937& rng, int count) -> std::vector<FaceVectorT> {
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<FaceVectorT> points;
  points.reserve(count);
  for (int i = 0; i < count; ++i) {
    double chi = dist(rng);
    double tau = dist(rng);
    if (chi + tau > 1.0) {
      chi = 1.0 - chi;
      tau = 1.0 - tau;
    }
    points.emplace_back(chi, tau);
  }
  return points;
}
} // namespace

TEST_CASE("Reference face map" * doctest::test_suite("geometry")) {
  constexpr double Epsilon = 1e-14;
  std::mt19937 rng(20260928);
  const auto points = samplePointsOnReferenceFace(rng, 200);

  SUBCASE("Lands on the face it names") {
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      for (const auto orientation : AllOrientations) {
        const ReferenceFaceMap map(side, orientation);
        for (const auto& point : points) {
          const auto cell = map.faceToCell(point);
          // the reference tetrahedron's faces are xi=0, eta=0, zeta=0 and xi+eta+zeta=1
          const std::array<double, Cell::NumFaces> coordinate{
              cell(2), cell(1), cell(0), cell(0) + cell(1) + cell(2) - 1.0};
          const std::array<std::size_t, Cell::NumFaces> zeroFor{0, 1, 2, 3};
          REQUIRE(coordinate[zeroFor[side]] == AbsApprox(0.0).epsilon(Epsilon));
        }
      }
    }
  }

  SUBCASE("Round trips face to cell and back") {
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      for (const auto orientation : AllOrientations) {
        const ReferenceFaceMap map(side, orientation);
        for (const auto& point : points) {
          const auto back = map.cellToFace(map.faceToCell(point));
          REQUIRE((back - point).norm() == AbsApprox(0.0).epsilon(Epsilon));
        }
      }
    }
  }

  SUBCASE("Its vertices are the face vertices of the reference cell") {
    const auto cellVertices = referenceCellVertices();
    const auto cellTransform = AffineTransform(cellVertices);

    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const AffineFaceTransform face(cellTransform, ReferenceFaceMap(side));
      const auto corners = face.vertices();
      for (std::size_t vertex = 0; vertex < Face::NumVertices; ++vertex) {
        const auto expected = cellVertices[MeshTools::FACE2NODES[side][vertex]];
        REQUIRE((corners[vertex] - expected).norm() == AbsApprox(0.0).epsilon(Epsilon));
      }
    }
  }
}

TEST_CASE("Affine face transform" * doctest::test_suite("geometry")) {
  constexpr double Epsilon = 1e-12;
  std::mt19937 rng(20260929);
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  // a well-shaped tetrahedron, so that the geometric quantities below are not dominated by
  // conditioning
  const std::array<CellVectorT, Cell::NumVertices> vertices{CellVectorT(0.1, 0.2, -0.3),
                                                            CellVectorT(2.1, 0.0, 0.4),
                                                            CellVectorT(-0.2, 1.7, 0.1),
                                                            CellVectorT(0.3, 0.1, 2.2)};
  const auto cell = AffineTransform(vertices);

  SUBCASE("Agrees with the cell transform on the face") {
    const auto points = samplePointsOnReferenceFace(rng, 100);
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      for (const auto orientation : AllOrientations) {
        const AffineFaceTransform face(cell, ReferenceFaceMap(side, orientation));
        for (const auto& point : points) {
          const auto viaFace = face.refToSpace(point);
          const auto viaCell = cell.refToSpace(face.refToCell(point));
          REQUIRE((viaFace - viaCell).norm() == AbsApprox(0.0).epsilon(Epsilon));
        }
      }
    }
  }

  SUBCASE("Reproduces the mesh normal, area and center") {
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const AffineFaceTransform face(cell, ReferenceFaceMap(side));
      const auto corners = face.vertices();

      // MeshTools::normal spans the face by its first two edges, in that order
      const CellVectorT expectedNormal = (corners[1] - corners[0]).cross(corners[2] - corners[0]);
      const auto normal = face.normal(FaceVectorT(Face::ReferenceBarycenter.data()));
      REQUIRE((normal - expectedNormal).norm() == AbsApprox(0.0).epsilon(Epsilon));

      REQUIRE(face.area() == AbsApprox(0.5 * expectedNormal.norm()).epsilon(Epsilon));

      const CellVectorT expectedCenter = (corners[0] + corners[1] + corners[2]) / 3.0;
      REQUIRE((face.center() - expectedCenter).norm() == AbsApprox(0.0).epsilon(Epsilon));
    }
  }

  SUBCASE("Its normal points away from the cell") {
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const AffineFaceTransform face(cell, ReferenceFaceMap(side));
      const auto normal = face.normal(FaceVectorT(Face::ReferenceBarycenter.data()));
      const CellVectorT toBarycenter =
          cell.refToSpace(CellVectorT(Cell::ReferenceBarycenter.data())) - face.center();
      REQUIRE(normal.dot(toBarycenter) < 0.0);
    }
  }

  SUBCASE("Its face-aligned basis is orthogonal and right-handed") {
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const AffineFaceTransform face(cell, ReferenceFaceMap(side));
      const auto basis = face.faceAlignedBasis();
      REQUIRE(basis[0].dot(basis[1]) == AbsApprox(0.0).epsilon(Epsilon));
      REQUIRE(basis[0].dot(basis[2]) == AbsApprox(0.0).epsilon(Epsilon));
      REQUIRE(basis[1].dot(basis[2]) == AbsApprox(0.0).epsilon(Epsilon));
      REQUIRE((basis[0].cross(basis[1]) - basis[2]).norm() == AbsApprox(0.0).epsilon(Epsilon));
    }
  }

  SUBCASE("Both sides of a shared face describe the same surface") {
    // Two cells sharing a triangle enumerate it differently. Combining the neighbor's side with
    // the matching orientation has to put the same face coordinate at the same point in space --
    // this is what lets quadrature points be matched across an interface, and it is the invariant
    // a wrong orientation table breaks.
    //
    // Here both cells carry the shared triangle on side 0, but list its vertices as (0, 2, 1) and
    // (1, 2, 0) respectively, which Rotate2 reconciles.
    const std::array<CellVectorT, Cell::NumVertices> neighborVertices{
        vertices[1], vertices[0], vertices[2], CellVectorT(0.4, 0.5, -2.4)};
    const auto neighbor = AffineTransform(neighborVertices);

    const AffineFaceTransform plus(cell, ReferenceFaceMap(0));
    const AffineFaceTransform minus(neighbor, ReferenceFaceMap(0, FaceOrientation::Rotate2));

    const auto points = samplePointsOnReferenceFace(rng, 100);
    for (const auto& point : points) {
      REQUIRE((plus.refToSpace(point) - minus.refToSpace(point)).norm() ==
              AbsApprox(0.0).epsilon(Epsilon));
    }

    // The orientation reconciles the parametrization, not the normal: it reverses the sense in
    // which the neighbor sweeps the face, so the reoriented normal no longer points away from the
    // neighbor. Each cell's outward normal comes from its own unoriented map, and those two are
    // antiparallel.
    const auto center = FaceVectorT(Face::ReferenceBarycenter.data());
    const AffineFaceTransform minusOwn(neighbor, ReferenceFaceMap(0));
    REQUIRE((plus.normal(center) + minusOwn.normal(center)).norm() ==
            AbsApprox(0.0).epsilon(Epsilon));
  }
}

} // namespace seissol::unit_test
