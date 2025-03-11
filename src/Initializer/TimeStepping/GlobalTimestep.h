// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_GLOBALTIMESTEP_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_GLOBALTIMESTEP_H_
#include <string>
#include <vector>

#include "Initializer/ParameterDB.h"

namespace seissol::initializer {
struct GlobalTimestep {
  std::vector<double> cellTimeStepWidths;
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
