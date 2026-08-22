// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "DynamicRupture/Output/Geometry.h"
#include "DynamicRupture/Output/OutputAux.h"

#include <array>
#include <cmath>

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

TEST_CASE("convertMaskFromBoolToInt" * doctest::test_suite("dynamicrupture")) {
  SUBCASE("Mixed") {
    const std::array<bool, 5> mask = {true, false, true, true, false};
    auto intMask = convertMaskFromBoolToInt<5>(mask);
    CHECK(intMask[0] == 1);
    CHECK(intMask[1] == 0);
    CHECK(intMask[2] == 1);
    CHECK(intMask[3] == 1);
    CHECK(intMask[4] == 0);
  }
  SUBCASE("All true") {
    const std::array<bool, 3> mask = {true, true, true};
    auto intMask = convertMaskFromBoolToInt<3>(mask);
    for (int i = 0; i < 3; ++i) {
      CHECK(intMask[i] == 1);
    }
  }
  SUBCASE("All false") {
    const std::array<bool, 3> mask = {false, false, false};
    auto intMask = convertMaskFromBoolToInt<3>(mask);
    for (int i = 0; i < 3; ++i) {
      CHECK(intMask[i] == 0);
    }
  }
}

} // namespace seissol::unit_test
