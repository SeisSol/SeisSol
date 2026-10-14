// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_INITIALIZER_FACESAMPLING_T_H_
#define SEISSOL_TESTS_INITIALIZER_FACESAMPLING_T_H_

// The material a face needs is asked for at the points of the face, each pulled into the cell
// it is asked for, so that the two cells sharing a face see their own sides of it.

#include <doctest.h>

#include "Common/Constants.h"
#include "Equations/Datastructures.h"
#include "Initializer/ParameterDB.h"
#include "TestHelper.h"

#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("Material at the points of a face" * doctest::test_suite("initializer")) {
  using namespace seissol::initializer;
  using Vertices = std::array<std::array<double, Cell::Dim>, Cell::NumVertices>;

  // Two cells sharing their face 0 (the vertices 0, 1, 2) in the plane y = 1/2, across which the
  // material of the test file jumps: below, a density of 2500, above, 1000, and on the plane
  // itself the one of below.
  const std::vector<Vertices> vertices{
      Vertices{{{0.0, 0.5, 0.0}, {1.0, 0.5, 0.0}, {0.0, 0.5, 1.0}, {0.2, 0.0, 0.3}}},
      Vertices{{{0.0, 0.5, 0.0}, {0.0, 0.5, 1.0}, {1.0, 0.5, 0.0}, {0.2, 1.0, 0.3}}},
  };
  const std::vector<int> groups{1, 1};
  const auto cells = CellToVertexArray::fromVectors(vertices, groups);

  const auto points = FacePointGenerator::fluxPoints();
  const FacePointGenerator generator(cells, points);
  REQUIRE(generator.outputPerCell() == Cell::NumFaces * points.count);
  const auto query = generator.generate();
  REQUIRE(query.numPoints() == vertices.size() * generator.outputPerCell());

  // the vertex each face is opposite to
  constexpr std::array<std::size_t, Cell::NumFaces> Opposite{3, 2, 1, 0};

  // Where the matrices of a face put its points, which is where the flux reads the field: the
  // point (chi, tau) of the reference triangle is (tau, chi, 0) on face 0, (chi, 0, tau) on face 1,
  // (0, tau, chi) on face 2 and (1 - chi - tau, chi, tau) on face 3. The material has to be asked
  // for at the same points in the same order.
  const auto onFace =
      [](std::size_t side,
         const std::array<double, Face::Dim>& point) -> std::array<double, Cell::Dim> {
    switch (side) {
    case 0:
      return {point[1], point[0], 0.0};
    case 1:
      return {point[0], 0.0, point[1]};
    case 2:
      return {0.0, point[1], point[0]};
    default:
      return {1.0 - point[0] - point[1], point[0], point[1]};
    }
  };

  for (std::size_t cell = 0; cell < vertices.size(); ++cell) {
    Eigen::Matrix3d edges;
    for (std::size_t j = 0; j < Cell::Dim; ++j) {
      for (std::size_t d = 0; d < Cell::Dim; ++d) {
        edges(d, j) = vertices[cell][j + 1][d] - vertices[cell][0][d];
      }
    }
    const Eigen::Vector3d origin(vertices[cell][0].data());
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      for (std::size_t i = 0; i < points.count; ++i) {
        const auto index = cell * generator.outputPerCell() + side * points.count + i;
        const Eigen::Vector3d point(query.x(index, 0), query.x(index, 1), query.x(index, 2));
        const Eigen::Vector3d reference = edges.inverse() * (point - origin);
        const std::array<double, Cell::NumVertices> barycentric{
            1.0 - reference.sum(), reference(0), reference(1), reference(2)};
        // inside the cell, and on its face but for the pull
        for (const auto coordinate : barycentric) {
          CHECK(coordinate > 0.0);
        }
        CHECK(barycentric[Opposite[side]] == doctest::Approx(InwardPull / Cell::NumVertices));
        const auto expected = onFace(side, points.point(i));
        for (std::size_t d = 0; d < Cell::Dim; ++d) {
          CHECK(reference(d) == doctest::Approx((1 - InwardPull) * expected[d] +
                                                InwardPull * Cell::ReferenceBarycenter[d])
                                    .epsilon(1e-10)
                                    .scale(1.0));
        }
        CHECK(query.group(index) == groups[cell]);
      }
    }
  }

  // each side of the shared face reads its own material there
  std::vector<seissol::model::MaterialT> materials;
  MaterialParameterDB<seissol::model::MaterialT> parameterDB;
  parameterDB.setMaterialVector(&materials);
  parameterDB.evaluateModel(tpath("Testing/material.yaml"), generator);
  REQUIRE(materials.size() == query.numPoints());
  for (std::size_t i = 0; i < points.count; ++i) {
    CHECK(materials[i].rho == doctest::Approx(2500.0));
    CHECK(materials[generator.outputPerCell() + i].rho == doctest::Approx(1000.0));
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_INITIALIZER_FACESAMPLING_T_H_
