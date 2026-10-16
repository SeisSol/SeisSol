// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GlobalTimestep.h"

#include "Common/Constants.h"
#include "Equations/Datastructures.h"
#include "Geometry/CellGeometry.h"
#include "Initializer/Model/WeightAdjusted.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters//SeisSolParameters.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Model/OperatorLayout.h"
#include "Parallel/MPI.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <mpi.h>
#include <utils/logger.h>
#include <vector>

namespace {

double computeCellTimestep(const std::array<Eigen::Vector3d, 4>& vertices,
                           double pWaveVel,
                           double cfl,
                           double maximumAllowedTimeStep) {
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

  // Every point a cell asks for the material at: where the material varies inside a cell, the
  // samples the operator is formed at, the points of the faces the flux is integrated at, and the
  // points of the plastic strain where those are others -- so that the timestep bounds every
  // wave speed a cell sees.
  std::vector<std::shared_ptr<seissol::initializer::QueryGenerator>> queryGens;
  // the sets the operator reads the material at, the first ones; whether a cell is weight-adjusted
  // follows from those
  std::size_t operatorSets = 0;
  if (seissolParams.model.materialNodal) {
    queryGens.push_back(std::make_shared<seissol::initializer::NodalPointGenerator>(
        cellToVertex, seissol::initializer::NodalPointGenerator::materialPoints()));
    if constexpr (NodalFlux) {
      queryGens.push_back(std::make_shared<seissol::initializer::FacePointGenerator>(
          cellToVertex, seissol::initializer::FacePointGenerator::fluxPoints()));
    }
    operatorSets = queryGens.size();
    if (seissolParams.model.plasticity &&
        !seissol::initializer::NodalPointGenerator::samePoints(
            seissol::initializer::NodalPointGenerator::plasticityPoints(),
            seissol::initializer::NodalPointGenerator::materialPoints())) {
      queryGens.push_back(std::make_shared<seissol::initializer::NodalPointGenerator>(
          cellToVertex, seissol::initializer::NodalPointGenerator::plasticityMaterialPoints()));
    }
  } else {
    queryGens.push_back(seissol::initializer::getBestQueryGenerator(
        seissolParams.model.useCellHomogenizedMaterial, cellToVertex));
  }

  std::vector<std::vector<Material>> materials(queryGens.size());
  std::size_t pointsPerCell = 0;
  for (std::size_t set = 0; set < queryGens.size(); ++set) {
    seissol::initializer::MaterialParameterDB<Material> parameterDB;
    parameterDB.setMaterialVector(&materials[set]);
    parameterDB.evaluateModel(seissolParams.model.materialFileName, *queryGens[set]);
    pointsPerCell += queryGens[set]->outputPerCell();
  }

  GlobalTimestep timestep;
  timestep.cellTimeStepWidths.resize(cellToVertex.size);
  timestep.weightAdjusted.assign(cellToVertex.size, 0);
  const double weightAdjustedRatio = seissolParams.model.weightAdjustedRatio;
  std::size_t localWeightAdjusted = 0;
  double localMinWeightAdjustedFactor = 1.0;

  // for the diagnostic below: the same timestep, but bounded by the mean wave
  // speed of a cell instead of its maximum
  double localMinMeanTimestep = std::numeric_limits<double>::max();
  double localMaxContrast = 1.0;

  for (unsigned cell = 0; cell < cellToVertex.size; ++cell) {
    // the wave speed bounds the cell from above, the material timestep from below
    double pWaveVel = 0.0;
    double pWaveVelSum = 0.0;
    double materialMaxTimestep = std::numeric_limits<double>::max();
    ImpedanceSpread spread;
    for (std::size_t set = 0; set < queryGens.size(); ++set) {
      const auto perCell = queryGens[set]->outputPerCell();
      for (std::size_t point = 0; point < perCell; ++point) {
        const auto& material = materials[set][cell * perCell + point];
        const auto pointWaveVel = material.getMaxWaveSpeed();
        pWaveVel = std::max(pWaveVel, pointWaveVel);
        pWaveVelSum += pointWaveVel;
        materialMaxTimestep = std::min(materialMaxTimestep, material.maximumTimestep());
        if (set < operatorSets) {
          spread.add(material);
        }
      }
    }
    // a weight-adjusted cell sees its wave speeds raised by what the lift without material
    // asks for (ImpedanceSpread::timestepFactor)
    double waveSpeedFactor = 1.0;
    if (isWeightAdjusted(spread.ratio(), weightAdjustedRatio)) {
      timestep.weightAdjusted[cell] = 1;
      ++localWeightAdjusted;
      localMinWeightAdjustedFactor =
          std::min(localMinWeightAdjustedFactor, spread.timestepFactor());
      waveSpeedFactor = 1.0 / spread.timestepFactor();
    }
    const double pWaveVelMean = pWaveVelSum / static_cast<double>(pointsPerCell);

    const std::array<Eigen::Vector3d, 4> vertices = cellToVertex.elementCoordinates(cell);
    const auto cellMaxTimestep =
        std::min(materialMaxTimestep, seissolParams.timeStepping.maxTimestepWidth);
    // a curved cell is as narrow as it is where it is thinnest, which its vertices do not say;
    // that scales what the wave speed permits, not the bounds from above
    const auto thickness = std::min(
        1.0, seissol::geometry::relativeThickness(*cellToVertex.elementTransform(cell), vertices));
    const auto cellTimestep = [&](double waveSpeed) {
      return std::min(thickness * computeCellTimestep(vertices,
                                                      waveSpeed,
                                                      seissolParams.timeStepping.cfl,
                                                      std::numeric_limits<double>::max()),
                      cellMaxTimestep);
    };
    timestep.cellTimeStepWidths[cell] = cellTimestep(pWaveVel * waveSpeedFactor);

    localMaxContrast = std::max(localMaxContrast, pWaveVel / pWaveVelMean);
    localMinMeanTimestep = std::min(localMinMeanTimestep, cellTimestep(pWaveVelMean));
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

  if (weightAdjustedRatio > 0.0) {
    std::size_t globalWeightAdjusted = 0;
    std::size_t globalCells = 0;
    std::size_t localCells = cellToVertex.size;
    double globalMinWeightAdjustedFactor = 1.0;
    MPI_Allreduce(&localWeightAdjusted,
                  &globalWeightAdjusted,
                  1,
                  Mpi::castToMpiType<std::size_t>(),
                  MPI_SUM,
                  seissol::Mpi::mpi.comm());
    MPI_Allreduce(&localCells,
                  &globalCells,
                  1,
                  Mpi::castToMpiType<std::size_t>(),
                  MPI_SUM,
                  seissol::Mpi::mpi.comm());
    MPI_Allreduce(&localMinWeightAdjustedFactor,
                  &globalMinWeightAdjustedFactor,
                  1,
                  MPI_DOUBLE,
                  MPI_MIN,
                  seissol::Mpi::mpi.comm());
    logInfo() << globalWeightAdjusted << "of" << globalCells
              << "cells are weight-adjusted: their impedances spread by more than"
              << weightAdjustedRatio << "- their time step is down to"
              << globalMinWeightAdjustedFactor << "times the one their wave speeds permit.";
  }

  if (pointsPerCell > 1) {
    double globalMaxContrast = 0;
    double globalMinMeanTimestep = 0;
    MPI_Allreduce(
        &localMaxContrast, &globalMaxContrast, 1, MPI_DOUBLE, MPI_MAX, seissol::Mpi::mpi.comm());
    MPI_Allreduce(&localMinMeanTimestep,
                  &globalMinMeanTimestep,
                  1,
                  MPI_DOUBLE,
                  MPI_MIN,
                  seissol::Mpi::mpi.comm());
    logInfo() << "Material asked for at" << pointsPerCell
              << "points per cell. Largest ratio between the maximum and the mean wave speed"
                 " of a cell:"
              << globalMaxContrast << "- the smallest timestep is"
              << timestep.globalMinTimeStep / globalMinMeanTimestep
              << "times the one the mean wave speeds would permit.";
  }

  return timestep;
}
} // namespace seissol::initializer
