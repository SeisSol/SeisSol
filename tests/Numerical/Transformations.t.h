// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Numerical/Transformation.h"

#include <Eigen/Dense>
#include <random>

namespace seissol::unit_test {

TEST_CASE("Test tetrahedron global to reference") {
  // We do all tests in double precision
  constexpr real Epsilon = 10 * std::numeric_limits<double>::epsilon();

  // NOLINTNEXTLINE (-cert-dcl59-cpp)
  std::mt19937 rnggen(9);
  std::uniform_real_distribution<> rngdist(0.0, 1.0);

  const auto vertices = std::array<Eigen::Vector3d, 4>{
      {Eigen::Vector3d(rngdist(rnggen), rngdist(rnggen), rngdist(rnggen)),
       Eigen::Vector3d(rngdist(rnggen), rngdist(rnggen), rngdist(rnggen)),
       Eigen::Vector3d(rngdist(rnggen), rngdist(rnggen), rngdist(rnggen)),
       Eigen::Vector3d(rngdist(rnggen), rngdist(rnggen), rngdist(rnggen))}};
  const auto center = 0.25 * (vertices[0] + vertices[1] + vertices[2] + vertices[3]);

  const auto res = seissol::transformations::tetrahedronGlobalToReference(
      vertices[0], vertices[1], vertices[2], vertices[3], center);
  REQUIRE(res(0) == AbsApprox(0.25).epsilon(Epsilon));
  REQUIRE(res(1) == AbsApprox(0.25).epsilon(Epsilon));
  REQUIRE(res(2) == AbsApprox(0.25).epsilon(Epsilon));
}

} // namespace seissol::unit_test
