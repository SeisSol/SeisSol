// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/PointMapper.h"
#include "Geometry/MockReader.h"

#include <Eigen/Dense>
#include <random>

namespace seissol::unit_test {

TEST_CASE("Point mapper") {
  // We do all tests in double precision
  std::array<Eigen::Vector3d, 4> vertices{};

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
  std::array<short, 3> contained{0, 0, 0};
  std::array<std::size_t, 3> meshId{std::numeric_limits<std::size_t>::max(),
                                    std::numeric_limits<std::size_t>::max(),
                                    std::numeric_limits<std::size_t>::max()};
  seissol::initializer::findMeshIds(points.data(), mockReader, 3, contained.data(), meshId.data());

  std::array<short, 3> expectedContained = {0, 1, 1};
  std::array<std::size_t, 3> expectedMeshId = {std::numeric_limits<std::size_t>::max(), 0, 0};

  REQUIRE(contained == expectedContained);
  REQUIRE(meshId == expectedMeshId);
}

} // namespace seissol::unit_test
