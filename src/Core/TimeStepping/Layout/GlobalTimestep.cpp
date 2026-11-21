// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GlobalTimestep.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>

namespace seissol::initializer {

double computeCellTimestep(const std::array<Eigen::Vector3d, 4>& vertices,
                           double pWaveVel,
                           double cfl,
                           double maximumAllowedTimeStep,
                           std::size_t convergenceOrder) {
  // Compute insphere radius
  std::array<Eigen::Vector3d, 4> x = vertices;
  Eigen::Matrix4d a;
  a << x[0](0), x[0](1), x[0](2), 1.0, x[1](0), x[1](1), x[1](2), 1.0, x[2](0), x[2](1), x[2](2),
      1.0, x[3](0), x[3](1), x[3](2), 1.0;

  const double alpha = a.determinant();
  const double nabc = (x[1] - x[0]).cross(x[2] - x[0]).norm();
  const double nabd = (x[1] - x[0]).cross(x[3] - x[0]).norm();
  const double nacd = (x[2] - x[0]).cross(x[3] - x[0]).norm();
  const double nbcd = (x[2] - x[1]).cross(x[3] - x[1]).norm();
  const double insphere = std::fabs(alpha) / (nabc + nabd + nacd + nbcd);

  // Compute maximum timestep
  return std::fmin(maximumAllowedTimeStep,
                   cfl * 2.0 * insphere / (pWaveVel * (2 * convergenceOrder - 1)));
}

} // namespace seissol::initializer
