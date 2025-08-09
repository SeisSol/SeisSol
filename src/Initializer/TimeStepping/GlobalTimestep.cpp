// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GlobalTimestep.h"

#include <Common/ConfigHelper.h>
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <mpi.h>
#include <variant>
#include <vector>

#include "Equations/Datastructures.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters//SeisSolParameters.h"

#include "Parallel/MPI.h"

namespace {

double
    computeCellTimestep(const std::array<Eigen::Vector3d, 4>& vertices,
                        double pWaveVel,
                        double cfl,
                        double maximumAllowedTimeStep,
                        std::size_t convergenceOrder,
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
                   cfl * 2.0 * insphere / (pWaveVel * (2 * convergenceOrder - 1)));
}

} // namespace

namespace seissol::initializer {

GlobalTimestep
    computeTimesteps(const seissol::initializer::CellToVertexArray& cellToVertex,
                     const seissol::initializer::parameters::SeisSolParameters& seissolParams) {

  const auto materials = queryMaterials(seissolParams.model, cellToVertex);

  GlobalTimestep timestep;
  timestep.cellTimeStepWidths.resize(cellToVertex.size);

  for (std::size_t cell = 0; cell < cellToVertex.size; ++cell) {
    std::visit(
        [&](auto cfg) {
          using Cfg = decltype(cfg);
          using MaterialT = model::MaterialTT<Cfg>;
          const auto& material = std::get<MaterialT>(materials[cell]);
          const double pWaveVel = material.getMaxWaveSpeed();
          const std::array<Eigen::Vector3d, 4> vertices = cellToVertex.elementCoordinates(cell);
          const auto materialMaxTimestep = material.maximumTimestep();
          const auto cellMaxTimestep =
              std::min(materialMaxTimestep, seissolParams.timeStepping.maxTimestepWidth);
          timestep.cellTimeStepWidths[cell] = computeCellTimestep(vertices,
                                                                  pWaveVel,
                                                                  seissolParams.timeStepping.cfl,
                                                                  cellMaxTimestep,
                                                                  Cfg::ConvergenceOrder,
                                                                  seissolParams);
        },
        ConfigVariantList[0]);
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
                seissol::MPI::mpi.comm());
  MPI_Allreduce(&localMaxTimestep,
                &timestep.globalMaxTimeStep,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                seissol::MPI::mpi.comm());
  return timestep;
}
} // namespace seissol::initializer
