// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Datastructures.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/LTSTree.h"
#include <Common/ConfigHelper.h>
#include <Common/Constants.h>
#include <Common/Real.h>
#include <Config.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/MemoryManager.h>
#include <Initializer/Parameters/ModelParameters.h>
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Kernels/Common.h>
#include <Kernels/Precision.h>
#include <Memory/Tree/Colormap.h>
#include <Memory/Tree/Layer.h>
#include <Model/CommonDatastructures.h>
#include <Model/Plasticity.h>
#include <Modules/Modules.h>
#include <Monitoring/Instrumentation.h>
#include <Monitoring/Stopwatch.h>
#include <Parallel/Helper.h>
#include <Physics/InstantaneousTimeMirrorManager.h>
#include <Solver/Estimator.h>
#include <Solver/MultipleSimulations.h>
#include <array>
#include <cassert>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <utils/env.h>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <variant>
#include <vector>

#include "InitModel.h"
#include "SeisSol.h"

using namespace seissol::initializer;

namespace {

using Plasticity = seissol::model::Plasticity;

void initializeCellMaterial(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  const auto& meshReader = seissolInstance.meshReader();
  initializer::MemoryManager& memoryManager = seissolInstance.getMemoryManager();

  // unpack ghost layer (merely a re-ordering operation, since the CellToVertexArray right now
  // requires an vector there)
  std::vector<std::array<std::array<double, Cell::Dim>, Cell::NumVertices>> ghostVertices;
  std::vector<int> ghostGroups;
  std::vector<int> ghostConfigs;
  std::unordered_map<int, std::vector<unsigned>> ghostIdxMap;
  for (const auto& neighbor : meshReader.getGhostlayerMetadata()) {
    ghostIdxMap[neighbor.first].reserve(neighbor.second.size());
    for (const auto& metadata : neighbor.second) {
      ghostIdxMap[neighbor.first].push_back(ghostVertices.size());
      std::array<std::array<double, Cell::Dim>, Cell::NumVertices> vertices{};
      for (size_t i = 0; i < Cell::NumVertices; ++i) {
        for (size_t j = 0; j < Cell::Dim; ++j) {
          vertices[i][j] = metadata.vertices[i][j];
        }
      }
      ghostVertices.emplace_back(vertices);
      ghostGroups.push_back(metadata.group);
      ghostConfigs.push_back(metadata.configId);
    }
  }

  const auto ghostOffset = meshReader.getElements().size();

  const auto meshArray = seissol::initializer::CellToVertexArray::fromMeshReader(meshReader);

  const auto meshHaloArray = seissol::initializer::CellToVertexArray::join(
      {meshArray,
       seissol::initializer::CellToVertexArray::fromVectors(
           ghostVertices, ghostGroups, ghostConfigs)});

  // material retrieval for ghost+copy+interior layers
  auto materials = queryMaterials(seissolParams.model, meshHaloArray);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (size_t i = 0; i < materials.size(); ++i) {
    auto& cellMat = materials[i];
    std::visit([&](auto& mat) { mat.initialize(seissolParams.model); }, cellMat);
  }

  // plasticity (if needed)
  std::vector<std::vector<Plasticity>> plasticityDB;
  if (seissolParams.model.plasticity) {
    // plasticity information is only needed on all interior+copy cells.
    const auto queryGen = seissol::initializer::getBestQueryGenerator<Plasticity>(
        seissolParams.model.plasticity, seissolParams.model.useCellHomogenizedMaterial, meshArray);

    const auto simcount = 1; // TODO (again)

    plasticityDB.resize(simcount);

    for (size_t i = 0; i < simcount; i++) {
      plasticityDB[i] = queryDB<Plasticity>(
          queryGen, seissolParams.model.plasticityFileNames[i], meshReader.getElements().size());
    }
  }

  logDebug() << "Setting cell materials in the storage (for interior and copy layers).";
  const auto& elements = meshReader.getElements();

  for (auto& layer : memoryManager.getLtsStorage().leaves()) {
    auto* cellInformation = layer.var<LTS::CellInformation>();
    auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();

    if (layer.getIdentifier().halo == HaloType::Ghost) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        const auto& localSecondaryInformation = secondaryInformation[cell];

        auto neighborRank = localSecondaryInformation.rank;
        auto neighborRankIdx = localSecondaryInformation.meshId;
        auto materialGhostIdx = ghostIdxMap.at(neighborRank)[neighborRankIdx];

        layer.wrap([&](auto cfg) {
          // NOLINTNEXTLINE
          auto& materialData = layer.var<LTS::MaterialData>(cfg)[cell];

          using Cfg = decltype(cfg);
          using MaterialT = model::MaterialTT<Cfg>;

          const auto& localMaterial =
              std::get<MaterialT>(materials[ghostOffset + materialGhostIdx]);

          initAssign(materialData, localMaterial);
        });
      }
    } else {
      auto* materialArray = layer.var<LTS::Material>();

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        // set the materials for the cell volume and its faces
        const auto& localSecondaryInformation = secondaryInformation[cell];
        const auto meshId = localSecondaryInformation.meshId;
        auto& material = materialArray[cell];
        const auto& localCellInformation = cellInformation[cell];

        layer.wrap([&](auto cfg) {
          // NOLINTNEXTLINE
          auto& materialData = layer.var<LTS::MaterialData>(cfg)[cell];

          using Cfg = decltype(cfg);
          using MaterialT = model::MaterialTT<Cfg>;

          const auto& localMaterial = std::get<MaterialT>(materials[meshId]);

          initAssign(materialData, localMaterial);

          material.local = &materialData;
        });

        for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
          if (isInternalFaceType(localCellInformation.faceTypes[side])) {
            // use the neighbor face material info in case that we are not at a boundary
            const auto& globalNeighborIndex = localSecondaryInformation.faceNeighbors[side];

            memoryManager.getLtsStorage().lookupWrap<LTS::MaterialData>(
                globalNeighborIndex,
                [&](auto& materialNeighbor) { material.neighbor[side] = &materialNeighbor; });
          } else {
            // otherwise, use the material from the own cell
            material.neighbor[side] = material.local;
          }
        }

        // if enabled, set up the plasticity as well
        if (seissolParams.model.plasticity) {
          const auto& localPlasticity = plasticityDB[meshId];

          layer.wrap([&](auto cfg) {
            using Cfg = decltype(cfg);

            auto& plasticity = layer.var<LTS::Plasticity>(cfg)[cell];
            assert(plasticityDB.size() == seissol::multisim::NumSimulations<Cfg> &&
                   "Plasticity database size mismatch with number of simulations");
            std::array<Plasticity, seissol::multisim::NumSimulations<Cfg>> localPlasticity;
            for (size_t i = 0; i < seissol::multisim::NumSimulations<Cfg>; ++i) {
              localPlasticity[i] = plasticityDB[i][meshId];
            }

            initAssign(plasticity,
                       seissol::model::PlasticityData<Cfg>(localPlasticity, material.local));
          });
        }
      }
    }
  }
}

void hostDeviceCoexecution(seissol::SeisSol& seissolInstance) {
  if constexpr (isDeviceOn()) {
    const auto hdswitch = seissolInstance.env().get<std::string>("DEVICE_HOST_SWITCH", "none");
    bool hdenabled = false;
    if (hdswitch == "none") {
      hdenabled = false;
      logInfo() << "No host-device switching. Everything runs on the GPU.";
    } else if (hdswitch == "auto") {
      hdenabled = true;
      logInfo() << "Automatic host-device switchpoint detection.";
      const auto hdswitchInt = solver::hostDeviceSwitch();
      seissolInstance.setExecutionPlaceCutoff(hdswitchInt);
    } else {
      hdenabled = true;
      const auto hdswitchInt = utils::StringUtils::parse<int>(hdswitch);
      logInfo() << "Manual host-device cutoff set to" << hdswitchInt << ".";
      seissolInstance.setExecutionPlaceCutoff(hdswitchInt);
    }

    const bool usmDefault = useUSM();

    if (!usmDefault && hdenabled) {
      logWarning() << "Using the host-device execution on non-USM systems is not fully supported "
                      "yet. Expect incorrect results.";
    }
  }
}

} // namespace

void seissol::initializer::initprocedure::initModel(seissol::SeisSol& seissolInstance) {
  SCOREP_USER_REGION("init_model", SCOREP_USER_REGION_TYPE_FUNCTION);

  logInfo() << "Begin init model.";

  // Call the pre mesh initialization hook
  seissol::Modules::callHook<ModuleHook::PreModel>();

  seissol::Stopwatch watch;
  watch.start();

  // these four methods need to be called in this order.
  logInfo() << "Model infos:";
  for (std::size_t i = 0; i < std::variant_size_v<ConfigVariant>; ++i) {
    logInfo() << "Model" << i;
    std::visit(
        [&](auto cfg) {
          using Cfg = decltype(cfg);
          logInfo() << "Material:" << model::MaterialTT<Cfg>::Text.c_str();
          logInfo() << "Order:" << Cfg::ConvergenceOrder;
          logInfo() << "Precision:"
                    << (Cfg::Precision == RealType::F32 ? "single (f32)" : "double (f64)");
          logInfo() << "Number of simulations: " << Cfg::NumSimulations;
        },
        ConfigVariantList[i]);
  }

  logInfo() << "Other model settings:";
  logInfo() << "Plasticity:"
            << (seissolInstance.getSeisSolParameters().model.plasticity ? "on" : "off");
  logInfo() << "Flux:"
            << parameters::fluxToString(seissolInstance.getSeisSolParameters().model.flux).c_str();
  logInfo() << "Flux near fault:"
            << parameters::fluxToString(seissolInstance.getSeisSolParameters().model.fluxNearFault)
                   .c_str();

  // init cell materials (needs LTS, to place the material in; this part was translated from
  // FORTRAN)
  logInfo() << "Initialize cell material parameters.";
  initializeCellMaterial(seissolInstance);

  if constexpr (isDeviceOn()) {
    logInfo() << "Determine Host-Device switchpoint";
    hostDeviceCoexecution(seissolInstance);
  }

  watch.pause();
  watch.printTime("Model initialized in:");

  // Call the post mesh initialization hook
  seissol::Modules::callHook<ModuleHook::PostModel>();

  logInfo() << "End init model.";
}
