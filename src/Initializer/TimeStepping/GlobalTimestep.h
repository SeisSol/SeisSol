// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_GLOBALTIMESTEP_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_GLOBALTIMESTEP_H_
#include <string>
#include <vector>

#include "Initializer/ParameterDB.h"

namespace seissol::initializer {
struct GlobalTimestep {
  std::vector<double> cellTimeStepWidths;
  double globalMinTimeStep;
  double globalMaxTimeStep;
};

namespace parameters {
struct SeisSolParameters;
} // namespace parameters

GlobalTimestep
    computeTimesteps(double cfl,
                     double maximumAllowedTimeStep,
                     const std::string& velocityModel,
                     const seissol::initializer::CellToVertexArray& cellToVertex,
                     const seissol::initializer::parameters::SeisSolParameters& seissolParams);
} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_GLOBALTIMESTEP_H_
