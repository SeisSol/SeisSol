#ifndef GLOBAL_TIMESTEP_HPP_
#define GLOBAL_TIMESTEP_HPP_
#include <vector>
#include <string>

#include "Initializer/ParameterDB.h"

namespace seissol::initializers {
struct GlobalTimestep {
  std::vector<double> cellTimeStepWidths;
  double globalMinTimeStep;
  double globalMaxTimeStep;
};

namespace parameters {
struct SeisSolParameters;
}

GlobalTimestep
    computeTimesteps(double cfl,
                     double maximumAllowedTimeStep,
                     const std::string& velocityModel,
                     const seissol::initializers::CellToVertexArray& cellToVertex,
                     const seissol::initializers::parameters::SeisSolParameters& seissolParams);
} // namespace seissol::initializers

#endif
