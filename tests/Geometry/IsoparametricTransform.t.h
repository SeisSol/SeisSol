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
#include "Geometry/IsoparametricTransform.h"
#include "TestHelper.h"

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <random>
#include <vector>

namespace seissol::unit_test::isoparametric {

using seissol::geometry::AffineTransform;
using seissol::geometry::CellTransform;
using seissol::geometry::FaceOrientation;
using seissol::geometry::IsoparametricFaceTransform;
using seissol::geometry::IsoparametricTransform;
using seissol::geometry::ReferenceFaceMap;
using VectorT = CellTransform::VectorEigenT;
using MatrixT = CellTransform::MatrixEigenT;
using FaceVectorT = seissol::geometry::FaceTransform::FaceVectorT;

/// a well-shaped tetrahedron, so that nothing below is dominated by conditioning
inline auto wellShaped(const VectorT& shift = VectorT::Zero())
    -> std::array<VectorT, Cell::NumVertices> {
  return {shift + VectorT(0.1, 0.2, -0.3),
          shift + VectorT(2.1, 0.0, 0.4),
          shift + VectorT(-0.2, 1.7, 0.1),
          shift + VectorT(0.3, 0.1, 2.2)};
}

/// the straight midpoints, moved by `amount` times a random displacement
inline auto bentMidpoints(const std::array<VectorT, Cell::NumVertices>& vertices,
                          std::mt19937& rng,
                          double amount) -> std::array<VectorT, 6> {
  constexpr std::array<std::array<std::size_t, 2>, 6> Edges{
      {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
  std::normal_distribution<double> gauss(0.0, 1.0);
  std::array<VectorT, 6> midpoints{};
  for (std::size_t edge = 0; edge < Edges.size(); ++edge) {
    midpoints[edge] = 0.5 * (vertices[Edges[edge][0]] + vertices[Edges[edge][1]]) +
                      amount * VectorT(gauss(rng), gauss(rng), gauss(rng));
  }
  return midpoints;
}

inline auto samplePointsInReferenceCell(std::mt19937& rng, int count) -> std::vector<VectorT> {
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<VectorT> points;
  while (points.size() < static_cast<std::size_t>(count)) {
    const VectorT point(dist(rng), dist(rng), dist(rng));
    if (point.sum() <= 1.0) {
      points.push_back(point);
    }
  }
  return points;
}

/// the cofactor matrix det(J) J^-1, whose rows a divergence is taken of
inline auto cofactor(const CellTransform& transform, const VectorT& point) -> MatrixT {
  const MatrixT jacobian = transform.refToSpaceJacobian(point);
  return jacobian.determinant() * jacobian.inverse();
}

TEST_CASE("Isoparametric cell transform" * doctest::test_suite("geometry")) {
  constexpr double Epsilon = 1e-12;
  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937 rng(20261006);
  const auto vertices = wellShaped();
  const auto points = samplePointsInReferenceCell(rng, 50);

  SUBCASE("The lattice holds the vertices first and as many nodes as the order asks for") {
    for (std::size_t order = 1; order <= 4; ++order) {
      const auto nodes = IsoparametricTransform::latticeNodes(order);
      REQUIRE(nodes.size() == (order + 1) * (order + 2) * (order + 3) / 6);
      REQUIRE((nodes[0] - VectorT(0, 0, 0)).norm() == AbsApprox(0.0).epsilon(Epsilon));
      REQUIRE((nodes[1] - VectorT(1, 0, 0)).norm() == AbsApprox(0.0).epsilon(Epsilon));
      REQUIRE((nodes[2] - VectorT(0, 1, 0)).norm() == AbsApprox(0.0).epsilon(Epsilon));
      REQUIRE((nodes[3] - VectorT(0, 0, 1)).norm() == AbsApprox(0.0).epsilon(Epsilon));
    }
  }

  SUBCASE("A straight-sided cell is the affine transform, whatever its order") {
    const auto affine = AffineTransform(vertices);
    for (std::size_t order = 1; order <= 3; ++order) {
      const auto transform = IsoparametricTransform::fromVertices(vertices, order);
      for (const auto& point : points) {
        REQUIRE((transform.refToSpace(point) - affine.refToSpace(point)).norm() ==
                AbsApprox(0.0).epsilon(Epsilon));
        REQUIRE((transform.refToSpaceJacobian(point) - affine.refToSpaceJacobian(point)).norm() ==
                AbsApprox(0.0).epsilon(Epsilon));
      }
    }
  }

  const auto curved =
      IsoparametricTransform::fromEdgeMidpoints(vertices, bentMidpoints(vertices, rng, 0.1));

  SUBCASE("A curved cell passes through its nodes") {
    const auto reference = IsoparametricTransform::latticeNodes(2);
    for (std::size_t i = 0; i < reference.size(); ++i) {
      REQUIRE((curved.refToSpace(reference[i]) - curved.nodes()[i]).norm() ==
              AbsApprox(0.0).epsilon(Epsilon));
    }
  }

  SUBCASE("The Jacobian is the derivative of the map") {
    constexpr double Step = 1e-6;
    for (const auto& point : points) {
      const MatrixT jacobian = curved.refToSpaceJacobian(point);
      for (std::size_t d = 0; d < Cell::Dim; ++d) {
        VectorT step = VectorT::Zero();
        step(d) = Step;
        const VectorT difference =
            (curved.refToSpace(VectorT(point + step)) - curved.refToSpace(VectorT(point - step))) /
            (2 * Step);
        REQUIRE((difference - jacobian.col(d)).norm() == AbsApprox(0.0).epsilon(1e-8));
      }
    }
  }

  SUBCASE("Its Jacobian varies inside the cell") {
    REQUIRE((curved.refToSpaceJacobian(VectorT(0.1, 0.1, 0.1)) -
             curved.refToSpaceJacobian(VectorT(0.6, 0.1, 0.1)))
                .norm() > 1e-3);
  }

  SUBCASE("The rows of the cofactor matrix are free of divergence") {
    // The metric identities: sum_e d/dxi_e (det(J) J^-1)_{ed} = 0. They are what lets the
    // divergence of a flux be written with the metric inside or outside the derivative, and a
    // discretisation that loses them does not keep a constant state constant.
    constexpr double Step = 1e-5;
    for (const auto& point : points) {
      VectorT divergence = VectorT::Zero();
      for (std::size_t e = 0; e < Cell::Dim; ++e) {
        VectorT step = VectorT::Zero();
        step(e) = Step;
        const MatrixT difference =
            (cofactor(curved, VectorT(point + step)) - cofactor(curved, VectorT(point - step))) /
            (2 * Step);
        divergence += difference.row(e).transpose();
      }
      REQUIRE(divergence.norm() == AbsApprox(0.0).epsilon(1e-7));
    }
  }

  SUBCASE("The Newton fallback inverts it") {
    for (const auto& point : points) {
      const auto back = curved.spaceToRef(curved.refToSpace(point));
      REQUIRE((back - point).norm() == AbsApprox(0.0).epsilon(1e-10));
    }
  }

  SUBCASE("Stays accurate for cells far away from the origin") {
    const auto shift = VectorT(5e5, 5e5, 0.0);
    const auto far = wellShaped(shift);
    std::array<VectorT, Cell::NumVertices> scaled{};
    for (std::size_t i = 0; i < Cell::NumVertices; ++i) {
      scaled[i] = shift + 100.0 * (far[i] - shift);
    }
    const auto transform =
        IsoparametricTransform::fromEdgeMidpoints(scaled, bentMidpoints(scaled, rng, 10.0));
    for (const auto& point : points) {
      const auto back = transform.spaceToRef(transform.refToSpace(point));
      REQUIRE((back - point).norm() == AbsApprox(0.0).epsilon(1e-9));
    }
  }
}

TEST_CASE("Isoparametric face transform" * doctest::test_suite("geometry")) {
  constexpr double Epsilon = 1e-12;
  // NOLINTNEXTLINE(bugprone-random-generator-seed,cert-msc32-c,cert-msc51-cpp)
  std::mt19937 rng(20261007);
  const auto vertices = wellShaped();
  const auto midpoints = bentMidpoints(vertices, rng, 0.1);
  const auto cell = IsoparametricTransform::fromEdgeMidpoints(vertices, midpoints);

  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::vector<FaceVectorT> points;
  while (points.size() < 50) {
    const FaceVectorT point(dist(rng), dist(rng));
    if (point.sum() <= 1.0) {
      points.push_back(point);
    }
  }

  SUBCASE("Is the cell's map on the face") {
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const IsoparametricFaceTransform face(cell, ReferenceFaceMap(side));
      for (const auto& point : points) {
        REQUIRE((face.refToSpace(point) - cell.refToSpace(face.refToCell(point))).norm() ==
                AbsApprox(0.0).epsilon(Epsilon));
      }
    }
  }

  SUBCASE("Its normal is the cofactor matrix applied to the reference normal") {
    // Nanson's relation, pointwise: (J a) x (J b) = det(J) J^-T (a x b). This is what makes the
    // surface Jacobian and the normal of a face agree with the metric of the cell at the face.
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const ReferenceFaceMap map(side);
      const IsoparametricFaceTransform face(cell, map);
      const VectorT referenceNormal =
          VectorT(map.embedding().col(0)).cross(VectorT(map.embedding().col(1)));
      for (const auto& point : points) {
        const MatrixT jacobian = cell.refToSpaceJacobian(face.refToCell(point));
        const VectorT expected =
            jacobian.determinant() * jacobian.inverse().transpose() * referenceNormal;
        REQUIRE((face.normal(point) - expected).norm() == AbsApprox(0.0).epsilon(1e-11));
      }
    }
  }

  SUBCASE("Its face-aligned basis turns with the normal") {
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const IsoparametricFaceTransform face(cell, ReferenceFaceMap(side));
      for (const auto& point : points) {
        const auto basis = face.faceAlignedBasis(point);
        const VectorT normal = face.normal(point);
        REQUIRE((basis[0] - normal).norm() == AbsApprox(0.0).epsilon(Epsilon));
        REQUIRE(basis[0].dot(basis[1]) == AbsApprox(0.0).epsilon(Epsilon));
        REQUIRE(basis[0].dot(basis[2]) == AbsApprox(0.0).epsilon(Epsilon));
        REQUIRE(basis[1].dot(basis[2]) == AbsApprox(0.0).epsilon(Epsilon));
        REQUIRE((basis[0].cross(basis[1]) - basis[2]).norm() == AbsApprox(0.0).epsilon(Epsilon));
      }
    }
    // and it does turn: the normals at two points of a curved face are not parallel
    const IsoparametricFaceTransform face(cell, ReferenceFaceMap(0));
    const VectorT first = face.normal(FaceVectorT(0.1, 0.1)).normalized();
    const VectorT second = face.normal(FaceVectorT(0.8, 0.1)).normalized();
    REQUIRE(first.cross(second).norm() > 1e-3);
  }

  SUBCASE("Both sides of a shared curved face describe the same surface") {
    // as in the affine case: the neighbour lists the shared triangle as (1, 0, 2), and Rotate2
    // reconciles the two parametrizations. The midpoints on the shared edges are the same points,
    // which is all a conforming curved mesh has to guarantee.
    const std::array<VectorT, Cell::NumVertices> neighborVertices{
        vertices[1], vertices[0], vertices[2], VectorT(0.4, 0.5, -2.4)};
    auto neighborMidpoints = bentMidpoints(neighborVertices, rng, 0.1);
    // edges (0,1), (0,2) and (1,2) of the neighbour are the edges (0,1), (1,2) and (0,2) here
    neighborMidpoints[0] = midpoints[0];
    neighborMidpoints[1] = midpoints[3];
    neighborMidpoints[3] = midpoints[1];
    const auto neighbor =
        IsoparametricTransform::fromEdgeMidpoints(neighborVertices, neighborMidpoints);

    const IsoparametricFaceTransform plus(cell, ReferenceFaceMap(0));
    const IsoparametricFaceTransform minus(neighbor, ReferenceFaceMap(0, FaceOrientation::Rotate2));
    const IsoparametricFaceTransform minusOwn(neighbor, ReferenceFaceMap(0));
    for (const auto& point : points) {
      REQUIRE((plus.refToSpace(point) - minus.refToSpace(point)).norm() ==
              AbsApprox(0.0).epsilon(Epsilon));
    }
    // and the outward normals of the two cells are antiparallel at every point of the face
    for (const auto& point : points) {
      const VectorT spacePoint = plus.refToSpace(point);
      const auto neighborPoint = minusOwn.map().cellToFace(neighbor.spaceToRef(spacePoint));
      REQUIRE((plus.normal(point) + minusOwn.normal(neighborPoint)).norm() ==
              AbsApprox(0.0).epsilon(1e-9));
    }
  }
}

} // namespace seissol::unit_test::isoparametric
