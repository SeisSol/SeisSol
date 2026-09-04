// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Instance/Geometry/Deduplicate.h"

#include <cmath>
#include <cstddef>
#include <random>
#include <set>
#include <vector>

namespace seissol::unit_test {

namespace {
using seissol::io::instance::geometry::deduplicatePoints;

//! The cube [0,1]^3 split into six tetrahedra, written the way the output writes them: every cell
//! repeats its four corners.
std::vector<double> cubeCorners() {
  const double vertices[8][3] = {
      {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
  const int cells[6][4] = {
      {0, 1, 3, 7}, {0, 1, 7, 5}, {0, 5, 7, 4}, {0, 3, 2, 7}, {0, 6, 4, 7}, {0, 2, 6, 7}};
  std::vector<double> result;
  for (const auto& cell : cells) {
    for (const auto vertex : cell) {
      result.insert(result.end(), &vertices[vertex][0], &vertices[vertex][0] + 3);
    }
  }
  return result;
}

} // namespace

TEST_CASE("IO/Deduplicate: coinciding corners become one point" * doctest::test_suite("io")) {
  const auto corners = cubeCorners();
  REQUIRE(corners.size() == 6 * 4 * 3);

  const auto map = deduplicatePoints(corners);

  // the eight corners of the cube, not the twenty-four the cells wrote
  CHECK(map.pointCount() == 8);
  REQUIRE(map.indices.size() == 24);

  // every input point still lands on its own coordinates
  for (std::size_t point = 0; point < map.indices.size(); ++point) {
    const auto target = map.indices[point];
    REQUIRE(target < map.pointCount());
    for (std::size_t d = 0; d < 3; ++d) {
      CHECK(map.points[target * 3 + d] == doctest::Approx(corners[point * 3 + d]));
    }
  }

  // and every unique point is actually reached
  const std::set<std::size_t> reached(map.indices.begin(), map.indices.end());
  CHECK(reached.size() == map.pointCount());
}

TEST_CASE("IO/Deduplicate: points that differ are kept apart" * doctest::test_suite("io")) {
  std::vector<double> points;
  std::mt19937 generator(42);
  std::uniform_real_distribution<double> position(-1.0, 1.0);
  for (std::size_t point = 0; point < 500; ++point) {
    for (std::size_t d = 0; d < 3; ++d) {
      points.push_back(position(generator));
    }
  }
  // every point once more, so the answer is known
  const auto original = points;
  points.insert(points.end(), original.begin(), original.end());

  const auto map = deduplicatePoints(points);
  CHECK(map.pointCount() == 500);
  for (std::size_t point = 0; point < 500; ++point) {
    CHECK(map.indices[point] == map.indices[point + 500]);
  }
}

TEST_CASE("IO/Deduplicate: a coordinate reached by two routes still merges" *
          doctest::test_suite("io")) {
  // The midpoint of an edge is computed once from each of the two cells that share it, and the
  // two results need not be bit-identical. Snapping alone would separate them whenever they fall
  // on either side of a cell boundary of the grid.
  const double a = 0.1 + 0.2;
  const double b = 0.3;
  REQUIRE(a != b);

  const std::vector<double> points{a, 1.0, 2.0, b, 1.0, 2.0, 0.5, 1.0, 2.0};
  const auto map = deduplicatePoints(points);

  CHECK(map.pointCount() == 2);
  CHECK(map.indices[0] == map.indices[1]);
  CHECK(map.indices[2] != map.indices[0]);
}

TEST_CASE("IO/Deduplicate: an empty input is handled" * doctest::test_suite("io")) {
  const auto map = deduplicatePoints({});
  CHECK(map.pointCount() == 0);
  CHECK(map.indices.empty());
}

} // namespace seissol::unit_test
