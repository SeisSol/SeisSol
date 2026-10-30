// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_GLOBALTIMESTEP_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_GLOBALTIMESTEP_H_
#include "Initializer/ParameterDB.h"

#include <cstdint>
#include <string>
#include <vector>

namespace seissol::initializer {
struct GlobalTimestep {
  std::vector<double> cellTimeStepWidths;
  // whether a cell is weight-adjusted (see WeightAdjustable), one byte per cell, so that the
  // partitioning can carry it along with the time step it went into
  std::vector<std::uint8_t> weightAdjusted;
  // whether a cell is curved (see initializer::CurvedCell), the same way; a curved cell is
  // weight-adjusted
  std::vector<std::uint8_t> curved;
  double globalMinTimeStep{std::numeric_limits<double>::max()};
  double globalMaxTimeStep{std::numeric_limits<double>::min()};
};

namespace parameters {
struct SeisSolParameters;
} // namespace parameters

GlobalTimestep
    computeTimesteps(const seissol::initializer::CellToVertexArray& cellToVertex,
                     const seissol::initializer::parameters::SeisSolParameters& seissolParams);
} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_GLOBALTIMESTEP_H_
