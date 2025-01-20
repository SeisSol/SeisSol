// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Numerical/Transformation.h"
#include <Eigen/Dense>

namespace seissol::unit_test {

TEST_CASE("Test tetrahedron global to reference") {
  // We do all tests in double precision
  constexpr real Epsilon = 10 * std::numeric_limits<double>::epsilon();

  std::srand(9);
  const auto vertices =
      std::array<Eigen::Vector3d, 4>{{Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX),
                                      Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX),
                                      Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX),
                                      Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX)}};
  const auto center = 0.25 * (vertices[0] + vertices[1] + vertices[2] + vertices[3]);

  const auto res = seissol::transformations::tetrahedronGlobalToReference(
      vertices[0].data(), vertices[1].data(), vertices[2].data(), vertices[3].data(), center);
  REQUIRE(res(0) == AbsApprox(0.25).epsilon(Epsilon));
  REQUIRE(res(1) == AbsApprox(0.25).epsilon(Epsilon));
  REQUIRE(res(2) == AbsApprox(0.25).epsilon(Epsilon));
}

} // namespace seissol::unit_test
