
#ifndef TIMESTEP_HPP_
#define TIMESTEP_HPP_
#include <Geometry/MeshReader.h>
#include <vector>
#include <string>

#include "Initializer/ParameterDB.h"

namespace seissol::initializer {
struct GlobalTimestep {
  std::vector<double> cellTimeStepWidths;
  double globalMinTimeStep;
  double globalMaxTimeStep;
};

struct LocalClustering {
  std::vector<int> cluster;
};

GlobalTimestep computeTimesteps(double cfl,
                                double maximumAllowedTimeStep,
                                const seissol::initializers::CellToVertexArray& cellToVertex);

std::vector<int> clusterTimesteps(const GlobalTimestep& timestep, double rate, double wiggle);

int enforceMaximumDifference(int maxDiff,
                             std::vector<int>& localClustering,
                             const geometry::MeshReader& meshReader);
} // namespace seissol::initializer

#endif
