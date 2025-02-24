// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/PointMapper.h"
#include "tests/Geometry/MockReader.h"
#include <Eigen/Dense>

namespace seissol::unit_test {

TEST_CASE("Point mapper") {
  // We do all tests in double precision
  std::array<Eigen::Vector3d, 4> vertices;

  std::srand(321);
  vertices = {{Eigen::Vector3d(0.0, 0.0, 0.0),
               Eigen::Vector3d(1.0, 0.0, 0.0),
               Eigen::Vector3d(0.0, 1.0, 0.0),
               Eigen::Vector3d(0.0, 0.0, 1.0)}};

  const seissol::MockReader mockReader(vertices);

  const Eigen::Vector3d points[3] = {Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                     (double)std::rand() / RAND_MAX,
                                                     (double)std::rand() / RAND_MAX),
                                     Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                     (double)std::rand() / RAND_MAX,
                                                     (double)std::rand() / RAND_MAX),
                                     0.25 *
                                         (vertices[0] + vertices[1] + vertices[2] + vertices[3])};
  short contained[3] = {0, 0, 0};
  unsigned meshId[3] = {std::numeric_limits<unsigned>::max(),
                        std::numeric_limits<unsigned>::max(),
                        std::numeric_limits<unsigned>::max()};
  seissol::initializer::findMeshIds(points, mockReader, 3, contained, meshId);

  std::array<short, 3> expectedContained = {0, 0, 1};
  std::array<unsigned, 3> expectedMeshId = {
      std::numeric_limits<unsigned>::max(), std::numeric_limits<unsigned>::max(), 0};

  for (int i = 0; i < 3; i++) {
    REQUIRE(contained[i] == expectedContained[i]);
    REQUIRE(meshId[i] == expectedMeshId[i]);
  }
}

} // namespace seissol::unit_test
