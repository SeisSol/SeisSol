// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_GLOBALTIMESTEP_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_GLOBALTIMESTEP_H_
#include "Initializer/ParameterDB.h"

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace seissol::initializer {
struct GlobalTimestep {
  std::vector<double> cellTimeStepWidths;
  double globalMinTimeStep{std::numeric_limits<double>::max()};
  double globalMaxTimeStep{std::numeric_limits<double>::min()};
};

namespace parameters {
struct SeisSolParameters;
} // namespace parameters

/**
 * The time step of a cell with the @p vertices and the largest wave speed @p pWaveVel by the CFL
 * condition with @p cfl for the order @p convergenceOrder , at most @p maximumAllowedTimeStep .
 */
double computeCellTimestep(const std::array<Eigen::Vector3d, 4>& vertices,
                           double pWaveVel,
                           double cfl,
                           double maximumAllowedTimeStep,
                           std::size_t convergenceOrder);

GlobalTimestep
    computeTimesteps(const seissol::initializer::CellToVertexArray& cellToVertex,
                     const seissol::initializer::parameters::SeisSolParameters& seissolParams);
} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_GLOBALTIMESTEP_H_
