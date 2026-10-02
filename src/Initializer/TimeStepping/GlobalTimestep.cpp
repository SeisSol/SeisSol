// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GlobalTimestep.h"

#include "Common/ConfigDispatch.h"
#include "Equations/Datastructures.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters//SeisSolParameters.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Parallel/MPI.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <mpi.h>
#include <vector>

namespace {

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
  const double nabc = ((x[1] - x[0]).cross(x[2] - x[0])).norm();
  const double nabd = ((x[1] - x[0]).cross(x[3] - x[0])).norm();
  const double nacd = ((x[2] - x[0]).cross(x[3] - x[0])).norm();
  const double nbcd = ((x[2] - x[1]).cross(x[3] - x[1])).norm();
  const double insphere = std::fabs(alpha) / (nabc + nabd + nacd + nbcd);

  // Compute maximum timestep
  return std::fmin(maximumAllowedTimeStep,
                   cfl * 2.0 * insphere / (pWaveVel * (2 * convergenceOrder - 1)));
}

} // namespace

namespace seissol::initializer {

GlobalTimestep
    computeTimesteps(const seissol::initializer::CellToVertexArray& cellToVertex,
                     const seissol::initializer::parameters::SeisSolParameters& seissolParams) {
  GlobalTimestep timestep;
  timestep.cellTimeStepWidths.resize(cellToVertex.size);

  // every cell with the material and the order of the configuration of its mesh group; the
  // material file is queried for the cells of each configuration separately, since another
  // material need not be defined in their groups
  const auto& model = seissolParams.model;
  const auto configs = model.configs();
  for (const auto config : configs) {
    std::vector<std::size_t> cells;
    for (std::size_t cell = 0; cell < cellToVertex.size; ++cell) {
      if (configs.size() == 1 || model.configOfGroup(cellToVertex.elementGroups(cell)) == config) {
        cells.push_back(cell);
      }
    }
    if (cells.empty()) {
      continue;
    }
    const auto cellsOfConfig =
        configs.size() == 1 ? cellToVertex
                            : seissol::initializer::CellToVertexArray::subset(cellToVertex, cells);

    dispatchConfig(config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      using Material = seissol::model::MaterialOf<Cfg>;

      const auto queryGen = seissol::initializer::getBestQueryGenerator<Material>(
          model.useCellHomogenizedMaterial, cellsOfConfig, Cfg::ConvergenceOrder);
      std::vector<Material> materials(cellsOfConfig.size);
      seissol::initializer::MaterialParameterDB<Material> parameterDB;
      parameterDB.setMaterialVector(&materials);
      parameterDB.evaluateModel(model.materialFileName, *queryGen);

      for (std::size_t i = 0; i < cellsOfConfig.size; ++i) {
        const double pWaveVel = materials[i].getMaxWaveSpeed();
        const std::array<Eigen::Vector3d, 4> vertices = cellsOfConfig.elementCoordinates(i);
        const auto materialMaxTimestep = materials[i].maximumTimestep();
        const auto cellMaxTimestep =
            std::min(materialMaxTimestep, seissolParams.timeStepping.maxTimestepWidth);
        timestep.cellTimeStepWidths[cells[i]] = computeCellTimestep(vertices,
                                                                    pWaveVel,
                                                                    seissolParams.timeStepping.cfl,
                                                                    cellMaxTimestep,
                                                                    Cfg::ConvergenceOrder);
      }
    });
  }

  const auto minmaxCellPosition =
      std::minmax_element(timestep.cellTimeStepWidths.begin(), timestep.cellTimeStepWidths.end());

  double localMinTimestep = *minmaxCellPosition.first;
  double localMaxTimestep = *minmaxCellPosition.second;

  MPI_Allreduce(&localMinTimestep,
                &timestep.globalMinTimeStep,
                1,
                MPI_DOUBLE,
                MPI_MIN,
                seissol::Mpi::mpi.comm());
  MPI_Allreduce(&localMaxTimestep,
                &timestep.globalMaxTimeStep,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                seissol::Mpi::mpi.comm());
  return timestep;
}
} // namespace seissol::initializer
