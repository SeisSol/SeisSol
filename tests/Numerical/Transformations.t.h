// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Geometry/CellTransform.h"
#include "Kernels/Precision.h"

#include <Eigen/Dense>
#include <array>
#include <limits>
#include <random>

namespace seissol::unit_test {

TEST_CASE("Test mapping a cell barycenter to reference coordinates" *
          doctest::test_suite("numerical")) {
  // We do all tests in double precision
  constexpr real Epsilon = 10 * std::numeric_limits<double>::epsilon();

  // NOLINTNEXTLINE (-cert-dcl59-cpp)
  std::mt19937 rnggen(9);
  std::uniform_real_distribution<> rngdist(0.0, 1.0);

  const auto vertices =
      std::array<CoordinateT, 4>{CoordinateT{rngdist(rnggen), rngdist(rnggen), rngdist(rnggen)},
                                 CoordinateT{rngdist(rnggen), rngdist(rnggen), rngdist(rnggen)},
                                 CoordinateT{rngdist(rnggen), rngdist(rnggen), rngdist(rnggen)},
                                 CoordinateT{rngdist(rnggen), rngdist(rnggen), rngdist(rnggen)}};
  const auto transform = seissol::geometry::AffineTransform(vertices);

  Eigen::Vector3d center = Eigen::Vector3d::Zero();
  for (const auto& vertex : vertices) {
    center += 0.25 * Eigen::Vector3d(vertex.data());
  }

  const auto res = transform.spaceToRef(center);
  CHECK(res(0) == AbsApprox(0.25).epsilon(Epsilon));
  CHECK(res(1) == AbsApprox(0.25).epsilon(Epsilon));
  CHECK(res(2) == AbsApprox(0.25).epsilon(Epsilon));
}

} // namespace seissol::unit_test
