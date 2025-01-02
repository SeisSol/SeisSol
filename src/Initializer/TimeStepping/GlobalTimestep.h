#ifndef GLOBAL_TIMESTEP_HPP_
#define GLOBAL_TIMESTEP_HPP_
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
    computeTimesteps(double cfl,
                     double maximumAllowedTimeStep,
                     const std::string& velocityModel,
                     const seissol::initializer::CellToVertexArray& cellToVertex,
                     const seissol::initializer::parameters::SeisSolParameters& seissolParams);
} // namespace seissol::initializer

#endif
