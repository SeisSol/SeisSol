// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Geometry/MockReader.h"
#include "Initializer/PointMapper.h"

#include <Eigen/Dense>
#include <random>

namespace seissol::unit_test {

TEST_CASE("Point mapper" * doctest::test_suite("initializer")) {
  // We do all tests in double precision
  std::array<Eigen::Vector3d, 4> vertices;

  vertices = {{Eigen::Vector3d(0.0, 0.0, 0.0),
               Eigen::Vector3d(1.0, 0.0, 0.0),
               Eigen::Vector3d(0.0, 1.0, 0.0),
               Eigen::Vector3d(0.0, 0.0, 1.0)}};

  const seissol::MockReader mockReader(vertices);

  // NOLINTNEXTLINE (-cert-dcl59-cpp)
  std::mt19937 rnggen(321);
  std::uniform_real_distribution<> rngdist(0.0, 1.0);

  const std::array<Eigen::Vector3d, 3> points{
      Eigen::Vector3d(rngdist(rnggen), rngdist(rnggen), rngdist(rnggen)),
      Eigen::Vector3d(rngdist(rnggen), rngdist(rnggen), rngdist(rnggen)),
      0.25 * (vertices[0] + vertices[1] + vertices[2] + vertices[3])};
  std::array<std::size_t, 3> meshId{std::numeric_limits<std::size_t>::max(),
                                    std::numeric_limits<std::size_t>::max(),
                                    std::numeric_limits<std::size_t>::max()};
  const auto contained =
      seissol::initializer::findUniqueMeshIds(points.data(), mockReader, 3, meshId.data());

  std::vector<bool> expectedContained{false, true, true};
  std::array<std::size_t, 3> expectedMeshId = {std::numeric_limits<std::size_t>::max(), 0, 0};

  CHECK(contained == expectedContained);
  CHECK(meshId == expectedMeshId);
}

TEST_CASE("Point mapper: a point on a shared face goes to the cell with the smaller global id" *
          doctest::test_suite("initializer")) {
  // two cells meeting in the face z = 0, both positively oriented
  const std::array<std::array<double, 3>, 5> coordinates{
      {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}}};
  std::vector<Vertex> vertices(coordinates.size());
  for (std::size_t i = 0; i < coordinates.size(); ++i) {
    std::copy(coordinates[i].begin(), coordinates[i].end(), vertices[i].coords);
  }
  const std::array<std::array<std::size_t, 4>, 2> cells{{{0, 1, 2, 3}, {0, 2, 1, 4}}};

  const std::array<Eigen::Vector3d, 2> points{Eigen::Vector3d(0.2, 0.2, 0.0),
                                              Eigen::Vector3d(0.1, 0.1, 0.1)};

  // the same two cells, in either order and with either of them holding the smaller global id
  for (const bool swapped : {false, true}) {
    for (const bool upperFirst : {false, true}) {
      std::vector<Element> elements(2);
      for (std::size_t local = 0; local < 2; ++local) {
        const auto cell = upperFirst ? local : 1 - local;
        std::copy(cells[cell].begin(), cells[cell].end(), elements[local].vertices);
        elements[local].localId = static_cast<LocalElemId>(local);
        elements[local].globalId = (cell == 0) == swapped ? 3 : 7;
      }
      const seissol::MockReader mockReader(vertices, elements);

      std::array<std::size_t, 2> meshId{std::numeric_limits<std::size_t>::max(),
                                        std::numeric_limits<std::size_t>::max()};
      const auto contained =
          seissol::initializer::findUniqueMeshIds(points.data(), mockReader, 2, meshId.data());

      REQUIRE(contained == std::vector<bool>{true, true});
      // on the shared face, both cells hold the point, and the smaller global id decides
      CHECK(elements[meshId[0]].globalId == 3);
      // inside the upper cell, only that one holds it
      const auto upper = upperFirst ? 0U : 1U;
      CHECK(meshId[1] == upper);
    }
  }
}

} // namespace seissol::unit_test
