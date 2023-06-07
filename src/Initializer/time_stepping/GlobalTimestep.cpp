#include "GlobalTimestep.hpp"

#include <vector>
#include <array>
#include <functional>
#include <Eigen/Dense>

#include "Equations/datastructures.hpp"
#include "Initializer/ParameterDB.h"
#include "Initializer/InputParameters.hpp"

#include "SeisSol.h"

namespace seissol::initializer {
static double computeCellTimestep(const std::array<Eigen::Vector3d, 4>& vertices,
                                  double pWaveVel,
                                  double cfl,
                                  double maximumAllowedTimeStep) {
  // Compute insphere radius
  std::array<Eigen::Vector3d, 4> x = vertices;
  Eigen::Matrix4d A;
  A << x[0](0), x[0](1), x[0](2), 1.0, x[1](0), x[1](1), x[1](2), 1.0, x[2](0), x[2](1), x[2](2),
      1.0, x[3](0), x[3](1), x[3](2), 1.0;

  double alpha = A.determinant();
  double Nabc = ((x[1] - x[0]).cross(x[2] - x[0])).norm();
  double Nabd = ((x[1] - x[0]).cross(x[3] - x[0])).norm();
  double Nacd = ((x[2] - x[0]).cross(x[3] - x[0])).norm();
  double Nbcd = ((x[2] - x[1]).cross(x[3] - x[1])).norm();
  double insphere = std::fabs(alpha) / (Nabc + Nabd + Nacd + Nbcd);

  // Compute maximum timestep
  return std::fmin(maximumAllowedTimeStep,
                   cfl * 2.0 * insphere / (pWaveVel * (2 * ConvergenceOrder - 1)));
}

GlobalTimestep computeTimesteps(double cfl,
                                double maximumAllowedTimeStep,
                                const std::string& velocityModel,
                                const seissol::initializers::CellToVertexArray& cellToVertex) {
  using Material = seissol::model::Material_t;

  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();

  auto* queryGen = seissol::initializers::getBestQueryGenerator(
      seissol::initializer::parameters::isModelAnelastic(),
      seissolParams.model.plasticity,
      seissol::initializer::parameters::isModelAnisotropic(),
      seissol::initializer::parameters::isModelPoroelastic(),
      seissolParams.model.useCellHomogenizedMaterial,
      cellToVertex);
  std::vector<Material> materials(cellToVertex.size);
  seissol::initializers::MaterialParameterDB<Material> parameterDB;
  parameterDB.setMaterialVector(&materials);
  parameterDB.evaluateModel(velocityModel, queryGen);

  GlobalTimestep timestep;
  timestep.cellTimeStepWidths.resize(cellToVertex.size);

  for (unsigned cell = 0; cell < cellToVertex.size; ++cell) {
    double pWaveVel = materials[cell].getMaxWaveSpeed();
    std::array<Eigen::Vector3d, 4> vertices = cellToVertex.elementCoordinates(cell);
    timestep.cellTimeStepWidths[cell] =
        computeCellTimestep(vertices, pWaveVel, cfl, maximumAllowedTimeStep);
  }

  const auto minmaxCellPosition =
      std::minmax_element(timestep.cellTimeStepWidths.begin(), timestep.cellTimeStepWidths.end());

  double localMinTimestep = *minmaxCellPosition.first;
  double localMaxTimestep = *minmaxCellPosition.second;

#ifdef USE_MPI
  MPI_Allreduce(&localMinTimestep,
                &timestep.globalMinTimeStep,
                1,
                MPI_DOUBLE,
                MPI_MIN,
                seissol::MPI::mpi.comm());
  MPI_Allreduce(&localMaxTimestep,
                &timestep.globalMaxTimeStep,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                seissol::MPI::mpi.comm());
#else
  timestep.globalMinTimeStep = localMinTimestep;
  timestep.globalMaxTimeStep = localMaxTimestep;
#endif
  return timestep;
}
} // namespace seissol::initializer
