// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GlobalTimestep.h"

#include <Common/Constants.h>
#include <Eigen/Dense>
#include <Initializer/Parameters/ModelParameters.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <mpi.h>
#include <vector>

#include "Equations/Datastructures.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters//SeisSolParameters.h"

namespace {

double
    computeCellTimestep(const std::array<Eigen::Vector3d, 4>& vertices,
                        double pWaveVel,
                        double cfl,
                        double maximumAllowedTimeStep,
                        const seissol::initializer::parameters::SeisSolParameters& seissolParams) {
  // Compute insphere radius
  std::array<Eigen::Vector3d, 4> x = vertices;
  Eigen::Matrix4d a;
  a << x[0](0), x[0](1), x[0](2), 1.0, x[1](0), x[1](1), x[1](2), 1.0, x[2](0), x[2](1), x[2](2),
      1.0, x[3](0), x[3](1), x[3](2), 1.0;

  const double alpha = a.determinant();
  const double nabc = ((x[1] - x[0]).cross(x[2] - x[0])).norm();
  const double nabd = ((x[1] - x[0]).cross(x[3] - x[0])).norm();
  const double nacd = ((x[2] - x[0]).cross(x[3] - x[0])).norm();
  const double nbcd = ((x[2] - x[1]).cross(x[3] - x[1])).norm();
  const double insphere = std::fabs(alpha) / (nabc + nabd + nacd + nbcd);

  // Compute maximum timestep
  return std::fmin(maximumAllowedTimeStep,
                   cfl * 2.0 * insphere / (pWaveVel * (2 * seissol::ConvergenceOrder - 1)));
}

} // namespace

namespace seissol::initializer {

GlobalTimestep
    computeTimesteps(const seissol::initializer::CellToVertexArray& cellToVertex,
                     const seissol::initializer::parameters::SeisSolParameters& seissolParams) {
  using Material = seissol::model::MaterialT;

  const auto queryGen = seissol::initializer::getBestQueryGenerator(
      seissolParams.model.plasticity, seissolParams.model.useCellHomogenizedMaterial, cellToVertex);
  std::vector<Material> materials(cellToVertex.size);
  seissol::initializer::MaterialParameterDB<Material> parameterDB;
  parameterDB.setMaterialVector(&materials);
  parameterDB.evaluateModel(seissolParams.model.materialFileName, *queryGen);

  GlobalTimestep timestep;
  timestep.cellTimeStepWidths.resize(cellToVertex.size);

  for (unsigned cell = 0; cell < cellToVertex.size; ++cell) {
    const double pWaveVel = materials[cell].getMaxWaveSpeed();
    const std::array<Eigen::Vector3d, 4> vertices = cellToVertex.elementCoordinates(cell);
    timestep.cellTimeStepWidths[cell] =
        computeCellTimestep(vertices,
                            pWaveVel,
                            seissolParams.timeStepping.cfl,
                            seissolParams.timeStepping.maxTimestepWidth,
                            seissolParams);
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
