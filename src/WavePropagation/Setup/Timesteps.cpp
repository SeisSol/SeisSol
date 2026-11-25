// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Timesteps.h"

#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/RunConfigs.h"
#include "Core/Data/DataReader.h"
#include "Core/Data/SimulationFiles.h"
#include "Core/Geometry/CellToVertexArray.h"
#include "Core/Parallel/MPI.h"
#include "Core/TimeStepping/Layout/GlobalTimestep.h"
#include "Core/TimeStepping/LtsParameters.h"
#include "WavePropagation/Physics/Datastructures.h"
#include "WavePropagation/Setup/MaterialParameterDB.h"
#include "WavePropagation/WavePropagationParameters.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <mpi.h>
#include <string>
#include <vector>

namespace seissol::initializer {

GlobalTimestep computeTimesteps(const CellToVertexArray& cellToVertex,
                                const parameters::RunConfigs& runConfigs,
                                const parameters::WavePropagationParameters& wavePropagation,
                                const parameters::TimeSteppingParameters& timeStepping,
                                const std::vector<std::string>& materialFileNames,
                                const data::ModelLoader& models) {
  GlobalTimestep timestep;
  timestep.cellTimeStepWidths.resize(cellToVertex.size);

  // every cell with the material and the order of the configuration of its mesh group; the
  // material file is queried for the cells of each configuration separately, since another
  // material need not be defined in their groups. The copies of a cell, which hold the simulations
  // of the run, step together, and so a cell takes the smallest time step that the materials of
  // the simulations allow.
  const auto configs = runConfigs.configs();
  const parameters::SimulationFiles materialFiles(materialFileNames);
  std::ranges::fill(timestep.cellTimeStepWidths, std::numeric_limits<double>::infinity());
  for (const auto config : configs) {
    std::vector<std::size_t> cells;
    for (std::size_t cell = 0; cell < cellToVertex.size; ++cell) {
      if (configs.size() == 1 ||
          runConfigs.configOfGroup(cellToVertex.elementGroups(cell)) == config) {
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
          wavePropagation.useCellHomogenizedMaterial, cellsOfConfig, Cfg::ConvergenceOrder);
      for (const auto& materialFile : materialFiles.files) {
        std::vector<Material> materials(cellsOfConfig.size);
        seissol::initializer::MaterialParameterDB<Material> parameterDB;
        parameterDB.setMaterialVector(&materials);
        parameterDB.evaluate(*models(materialFile, {"x", "y", "z"}), *queryGen);

        for (std::size_t i = 0; i < cellsOfConfig.size; ++i) {
          const double pWaveVel = materials[i].getMaxWaveSpeed();
          const std::array<Eigen::Vector3d, 4> vertices = cellsOfConfig.elementCoordinates(i);
          const auto materialMaxTimestep = materials[i].maximumTimestep();
          const auto cellMaxTimestep = std::min(materialMaxTimestep, timeStepping.maxTimestepWidth);
          timestep.cellTimeStepWidths[cells[i]] = std::min(
              timestep.cellTimeStepWidths[cells[i]],
              computeCellTimestep(
                  vertices, pWaveVel, timeStepping.cfl, cellMaxTimestep, Cfg::ConvergenceOrder));
        }
      }
    });
  }

  const auto minmaxCellPosition = std::ranges::minmax_element(timestep.cellTimeStepWidths);

  const double localMinTimestep = *minmaxCellPosition.min;
  const double localMaxTimestep = *minmaxCellPosition.max;

  timestep.globalMinTimeStep = seissol::Mpi::mpi.allreduce(localMinTimestep, MPI_MIN);
  timestep.globalMaxTimeStep = seissol::Mpi::mpi.allreduce(localMaxTimestep, MPI_MAX);
  return timestep;
}

std::string ltsUnsupportedMaterial(const parameters::RunConfigs& runConfigs) {
  std::string material;
  for (const auto config : runConfigs.configs()) {
    dispatchConfig(config, [&](auto cfg) {
      using MaterialT = model::MaterialOf<decltype(cfg)>;
      if (!MaterialT::SupportsLTS && material.empty()) {
        material = MaterialT::Text;
      }
    });
  }
  return material;
}

} // namespace seissol::initializer
