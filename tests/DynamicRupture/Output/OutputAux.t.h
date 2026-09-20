// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Constants.h"
#include "DynamicRupture/Output/Geometry.h"
#include "DynamicRupture/Output/OutputAux.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::dr;

TEST_CASE("ExtVrtxCoords construction" * doctest::test_suite("dynamicrupture")) {
  SUBCASE("Default zero") {
    ExtVrtxCoords p;
    CHECK(p[0] == doctest::Approx(0.0));
    CHECK(p[1] == doctest::Approx(0.0));
    CHECK(p[2] == doctest::Approx(0.0));
  }
  SUBCASE("Initializer list") {
    ExtVrtxCoords p = {1.0, 2.0, 3.0};
    CHECK(p[0] == doctest::Approx(1.0));
    CHECK(p[1] == doctest::Approx(2.0));
    CHECK(p[2] == doctest::Approx(3.0));
    CHECK(ExtVrtxCoords::size() == 3);
  }
  SUBCASE("Eigen conversion") {
    ExtVrtxCoords p = {4.0, 5.0, 6.0};
    auto v = p.getAsEigen3LibVector();
    CHECK(v[0] == doctest::Approx(4.0));
    CHECK(v[1] == doctest::Approx(5.0));
    CHECK(v[2] == doctest::Approx(6.0));
  }
}

TEST_CASE("ExtTriangle construction" * doctest::test_suite("dynamicrupture")) {
  ExtTriangle tri(
      ExtVrtxCoords{0.0, 0.0, 0.0}, ExtVrtxCoords{1.0, 0.0, 0.0}, ExtVrtxCoords{0.0, 1.0, 0.0});
  CHECK(tri.point(0)[0] == doctest::Approx(0.0));
  CHECK(tri.point(1)[0] == doctest::Approx(1.0));
  CHECK(tri.point(2)[1] == doctest::Approx(1.0));
  CHECK(ExtTriangle::size() == 3);
}

TEST_CASE("getReferenceTriangle" * doctest::test_suite("dynamicrupture")) {
  for (int side = 0; side < 4; ++side) {
    auto tri = getReferenceTriangle(side);
    for (int v = 0; v < 3; ++v) {
      for (int d = 0; d < 3; ++d) {
        CHECK(tri.point(v)[d] >= -1e-15);
        CHECK(tri.point(v)[d] <= 1.0 + 1e-15);
      }
    }
  }
  SUBCASE("Side 0 in z=0 plane") {
    auto tri = getReferenceTriangle(0);
    CHECK(tri.point(0)[2] == doctest::Approx(0.0));
    CHECK(tri.point(1)[2] == doctest::Approx(0.0));
    CHECK(tri.point(2)[2] == doctest::Approx(0.0));
  }
  SUBCASE("Side 3 oblique") {
    auto tri = getReferenceTriangle(3);
    CHECK(tri.point(0)[0] == doctest::Approx(1.0));
    CHECK(tri.point(1)[1] == doctest::Approx(1.0));
    CHECK(tri.point(2)[2] == doctest::Approx(1.0));
  }
}

TEST_CASE("getMidPointTriangle" * doctest::test_suite("dynamicrupture")) {
  ExtTriangle tri(
      ExtVrtxCoords{0.0, 0.0, 0.0}, ExtVrtxCoords{3.0, 0.0, 0.0}, ExtVrtxCoords{0.0, 6.0, 0.0});
  auto mid = getMidPointTriangle(tri);
  CHECK(mid[0] == doctest::Approx(1.0));
  CHECK(mid[1] == doctest::Approx(2.0));
  CHECK(mid[2] == doctest::Approx(0.0));
}

TEST_CASE("getMidPoint" * doctest::test_suite("dynamicrupture")) {
  ExtVrtxCoords a = {2.0, 4.0, 6.0};
  ExtVrtxCoords b = {8.0, 10.0, 12.0};
  auto mid = getMidPoint(a, b);
  CHECK(mid[0] == doctest::Approx(5.0));
  CHECK(mid[1] == doctest::Approx(7.0));
  CHECK(mid[2] == doctest::Approx(9.0));
}

TEST_CASE("computeTriangleArea" * doctest::test_suite("dynamicrupture")) {
  SUBCASE("Unit right triangle") {
    ExtTriangle tri(
        ExtVrtxCoords{0.0, 0.0, 0.0}, ExtVrtxCoords{1.0, 0.0, 0.0}, ExtVrtxCoords{0.0, 1.0, 0.0});
    CHECK(computeTriangleArea(tri) == doctest::Approx(0.5));
  }
  SUBCASE("Equilateral triangle side=2") {
    ExtTriangle tri(ExtVrtxCoords{0.0, 0.0, 0.0},
                    ExtVrtxCoords{2.0, 0.0, 0.0},
                    ExtVrtxCoords{1.0, std::sqrt(3.0), 0.0});
    CHECK(computeTriangleArea(tri) == doctest::Approx(std::sqrt(3.0)));
  }
  SUBCASE("3D triangle") {
    ExtTriangle tri(
        ExtVrtxCoords{0.0, 0.0, 0.0}, ExtVrtxCoords{1.0, 0.0, 0.0}, ExtVrtxCoords{0.0, 0.0, 1.0});
    CHECK(computeTriangleArea(tri) == doctest::Approx(0.5));
  }
  SUBCASE("Degenerate") {
    ExtTriangle tri(
        ExtVrtxCoords{0.0, 0.0, 0.0}, ExtVrtxCoords{1.0, 0.0, 0.0}, ExtVrtxCoords{2.0, 0.0, 0.0});
    CHECK(computeTriangleArea(tri) == doctest::Approx(0.0).epsilon(1e-15));
  }
}

TEST_CASE("getDistanceFromPointToFace" * doctest::test_suite("dynamicrupture")) {
  ExtTriangle face(
      ExtVrtxCoords{0.0, 0.0, 0.0}, ExtVrtxCoords{1.0, 0.0, 0.0}, ExtVrtxCoords{0.0, 1.0, 0.0});
  VrtxCoords normal = {0.0, 0.0, 1.0};
  SUBCASE("Above") {
    ExtVrtxCoords pt = {0.5, 0.5, 3.0};
    CHECK(getDistanceFromPointToFace(pt, face, normal) == doctest::Approx(-3.0));
  }
  SUBCASE("On") {
    ExtVrtxCoords pt = {0.25, 0.25, 0.0};
    CHECK(getDistanceFromPointToFace(pt, face, normal) == doctest::Approx(0.0));
  }
  SUBCASE("Below") {
    ExtVrtxCoords pt = {0.5, 0.5, -2.0};
    CHECK(getDistanceFromPointToFace(pt, face, normal) == doctest::Approx(2.0));
  }
}

TEST_CASE("projectPointToFace" * doctest::test_suite("dynamicrupture")) {
  ExtTriangle face(
      ExtVrtxCoords{0.0, 0.0, 0.0}, ExtVrtxCoords{1.0, 0.0, 0.0}, ExtVrtxCoords{0.0, 1.0, 0.0});
  VrtxCoords normal = {0.0, 0.0, 1.0};
  ExtVrtxCoords pt = {0.3, 0.3, 5.0};
  projectPointToFace(pt, face, normal);
  CHECK(pt[0] == doctest::Approx(0.3));
  CHECK(pt[1] == doctest::Approx(0.3));
  CHECK(pt[2] == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("getNearestFacePoint" * doctest::test_suite("dynamicrupture")) {
  double facePoints[][2] = {{0.0, 0.0}, {1.0, 0.0}, {0.5, 0.5}, {0.0, 1.0}};
  SUBCASE("Exact match") {
    double target[2] = {1.0, 0.0};
    auto [idx, dist] = getNearestFacePoint(target, facePoints, 4);
    CHECK(idx == 1);
    CHECK(dist == doctest::Approx(0.0));
  }
  SUBCASE("Nearest to center") {
    double target[2] = {0.4, 0.4};
    auto [idx, dist] = getNearestFacePoint(target, facePoints, 4);
    CHECK(idx == 2);
  }
  SUBCASE("Nearest to origin") {
    double target[2] = {0.01, 0.01};
    auto [idx, dist] = getNearestFacePoint(target, facePoints, 4);
    CHECK(idx == 0);
  }
}

TEST_CASE("getClosestInternalStroudGp" * doctest::test_suite("dynamicrupture")) {
  SUBCASE("Interior stays") { CHECK(getClosestInternalStroudGp(15, 4) == 15); }
  SUBCASE("Edge moves inward") { CHECK(getClosestInternalStroudGp(3, 4) == 9); }
}

TEST_CASE("getElementVertexId ranges" * doctest::test_suite("dynamicrupture")) {
  for (int side = 0; side < 4; ++side) {
    std::array<bool, 4> used = {false, false, false, false};
    for (int fv = 0; fv < 3; ++fv) {
      int vid = getElementVertexId(side, fv);
      CHECK(vid >= 0);
      CHECK(vid <= 3);
      used[vid] = true;
    }
    int count = 0;
    for (const bool u : used) {
      if (u) {
        ++count;
      }
    }
    CHECK(count == 3);
  }
}

//! A fault output of @p cellCount cells, emitted the way the refiner emits it: the receivers of a
//! cell lie consecutively, point-major and simulation-minor.
inline ReceiverPoints makeFaultOutput(std::size_t cellCount,
                               std::size_t pointsPerCell,
                               std::size_t simulationCount,
                               const std::vector<int>& tags,
                               const std::vector<std::size_t>& elements,
                               const std::vector<int>& sides) {
  ReceiverPoints points(cellCount * pointsPerCell * simulationCount);
  for (std::size_t cell = 0; cell < cellCount; ++cell) {
    for (std::size_t point = 0; point < pointsPerCell; ++point) {
      for (std::size_t sim = 0; sim < simulationCount; ++sim) {
        auto& receiver = points[(cell * pointsPerCell + point) * simulationCount + sim];
        receiver.faultTag = tags[cell];
        receiver.elementGlobalIndex = elements[cell];
        receiver.localFaceSideId = sides[cell];
        receiver.simIndex = static_cast<int>(sim);
      }
    }
  }
  return points;
}

TEST_CASE("fault output cell properties" * doctest::test_suite("dynamicrupture")) {
  // Two faces of the same element. The tags repeat across faces and are chosen so that no tag
  // coincides with a face identifier: a cell property read from the wrong field cannot pass.
  const std::vector<int> tags{101, 101};
  const std::vector<std::size_t> elements{7, 7};
  const std::vector<int> sides{2, 3};

  SUBCASE("Order zero, one simulation") {
    const auto points = makeFaultOutput(2, 1, 1, tags, elements, sides);

    CHECK(faultTagOfCell(points, 0, 1, 1) == 101);
    CHECK(faultTagOfCell(points, 1, 1, 1) == 101);

    CHECK(globalFaceIdOfCell(points, 0, 1, 1) == 7 * Cell::NumFaces + 2);
    CHECK(globalFaceIdOfCell(points, 1, 1, 1) == 7 * Cell::NumFaces + 3);
  }

  SUBCASE("Higher order") {
    const std::size_t pointsPerCell = 3;
    const auto points = makeFaultOutput(2, pointsPerCell, 1, tags, elements, sides);

    CHECK(firstReceiverOfCell(1, pointsPerCell, 1) == pointsPerCell);
    CHECK(faultTagOfCell(points, 1, pointsPerCell, 1) == 101);
    CHECK(globalFaceIdOfCell(points, 1, pointsPerCell, 1) == 7 * Cell::NumFaces + 3);
  }

  SUBCASE("Fused simulations") {
    const std::size_t pointsPerCell = 3;
    const std::size_t simulationCount = 4;
    const auto points = makeFaultOutput(2, pointsPerCell, simulationCount, tags, elements, sides);

    CHECK(firstReceiverOfCell(1, pointsPerCell, simulationCount) ==
          pointsPerCell * simulationCount);
    CHECK(faultTagOfCell(points, 1, pointsPerCell, simulationCount) == 101);
    CHECK(globalFaceIdOfCell(points, 1, pointsPerCell, simulationCount) == 7 * Cell::NumFaces + 3);
  }

  SUBCASE("The tag is not the face identifier") {
    const auto points = makeFaultOutput(2, 1, 1, tags, elements, sides);

    // Two cells that share a tag but not a face: a field that told them apart would not be the
    // tag, and one that did not would not be the identifier.
    CHECK(faultTagOfCell(points, 0, 1, 1) == faultTagOfCell(points, 1, 1, 1));
    CHECK(globalFaceIdOfCell(points, 0, 1, 1) != globalFaceIdOfCell(points, 1, 1, 1));

    for (std::size_t cell = 0; cell < 2; ++cell) {
      CHECK(static_cast<std::size_t>(faultTagOfCell(points, cell, 1, 1)) !=
            globalFaceIdOfCell(points, cell, 1, 1));
    }
  }
}

} // namespace seissol::unit_test
