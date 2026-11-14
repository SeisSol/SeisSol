// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitSideConditions.h"

#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Initializer/InitialFieldProjection.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Parameters/SimulationFiles.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/LTS.h"
#include "Physics/InitialField.h"
#include "Physics/Scenario/Registry.h"
#include "SeisSol.h"
#include "SourceTerm/Manager.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer::initprocedure {

namespace {

/// The initial conditions of the cells of the configuration `config`, set up with the material
/// `materialData` of one of them; for the simulation `simulation` alone, if given.
std::vector<std::unique_ptr<physics::InitialField>>
    buildInitialConditionList(seissol::SeisSol& seissolInstance,
                              ConfigId config,
                              const CellMaterialData& materialData,
                              std::optional<std::size_t> simulation = {}) {
  const auto& parameters = seissolInstance.parameters();
  const auto type = parameters.initialization.type;

  const auto availability = physics::scenario::availability(type, config);
  if (!availability.available) {
    logError() << "The initial condition" << physics::scenario::name(type).data()
               << "cannot be used with material"
               << materialTypeName(configValue(config).materialType).data() << "--"
               << availability.reason.data() << ".";
  }

  return physics::scenario::build(
      type,
      physics::scenario::Input{.parameters = parameters,
                               .materialData = materialData,
                               .gravitation = seissolInstance.gravitationSetup(),
                               .config = config,
                               .simulation = simulation});
}

/// The initial conditions of the cells of the configuration `config`, set up with the material of
/// the mesh cell `cell`. Where the copies of the cells differ in their material, a simulation is
/// set up with the one of the copy of `cell` that holds it, and so on its own; otherwise, the
/// simulations are set up together, with the material of the first copy.
std::vector<std::unique_ptr<physics::InitialField>>
    buildInitialConditions(seissol::SeisSol& seissolInstance, ConfigId config, std::size_t cell) {
  auto& memoryManager = seissolInstance.memoryManager();
  auto& storage = memoryManager.ltsStorage();
  const auto& backmap = memoryManager.backmap();
  const auto& settings = memoryManager.simulationSettings();
  const auto& materialFiles = seissolInstance.parameters().model.materialFileNames;

  bool copiesDiffer = false;
  for (std::size_t copy = 1; copy < settings.copyCount(); ++copy) {
    copiesDiffer = copiesDiffer || materialFiles[copy * settings.fusedCount] != materialFiles[0];
  }
  if (!copiesDiffer) {
    return buildInitialConditionList(
        seissolInstance, config, storage.lookup<LTS::Material>(backmap.get(cell, 0)));
  }

  std::vector<std::unique_ptr<physics::InitialField>> fields;
  fields.reserve(settings.simulationCount);
  for (std::size_t simulation = 0; simulation < settings.simulationCount; ++simulation) {
    const auto copy = simulation / settings.fusedCount;
    auto field = buildInitialConditionList(seissolInstance,
                                           config,
                                           storage.lookup<LTS::Material>(backmap.get(cell, copy)),
                                           simulation);
    fields.push_back(std::move(field.at(0)));
  }
  return fields;
}

void initInitialCondition(seissol::SeisSol& seissolInstance) {
  const auto& initConditionParams = seissolInstance.parameters().initialization;
  auto& memoryManager = seissolInstance.memoryManager();

  if (initConditionParams.type == seissol::initializer::parameters::InitializationType::Easi) {
    const parameters::SimulationFiles iniFields(initConditionParams.filenames);
    for (const auto& file : iniFields.files) {
      logInfo() << "Loading the initial condition from the easi file" << file;
    }
    seissol::initializer::projectEasiInitialField(iniFields,
                                                  seissolInstance.meshReader(),
                                                  memoryManager.ltsStorage(),
                                                  initConditionParams.hasTime);
  } else {
    logInfo() << "Using initial condition"
              << physics::scenario::name(initConditionParams.type).data() << ".";

    // The configurations of one material set up the scenario alike, with the material of the
    // first cell of the first of them in the run that has cells on this rank: setting up a
    // scenario may pick among degenerate eigenvectors of the material, which then roundoff
    // decides, and configurations can average the material of a cell differently (e.g. with the
    // quadrature of another order).
    auto& storage = memoryManager.ltsStorage();
    const auto& backmap = memoryManager.backmap();
    const auto cellCount = seissolInstance.meshReader().getElements().size();
    std::vector<std::optional<std::size_t>> firstCell(builtConfigCount());
    for (std::size_t cell = 0; cell < cellCount; ++cell) {
      const auto config = storage.lookup<LTS::SecondaryInformation>(backmap.get(cell, 0)).configId;
      if (!firstCell[config].has_value()) {
        firstCell[config] = cell;
      }
    }
    const auto sameMaterial = [](ConfigId first, ConfigId second) {
      return configValue(first).materialType == configValue(second).materialType &&
             configValue(first).relaxationMechanisms == configValue(second).relaxationMechanisms;
    };
    const auto configs = seissolInstance.parameters().model.configs();
    for (const auto config : configs) {
      for (const auto reference : configs) {
        if (sameMaterial(reference, config) && firstCell[reference].has_value()) {
          memoryManager.setInitialConditions(
              config,
              buildInitialConditions(seissolInstance, config, firstCell[reference].value()));
          break;
        }
      }
    }

    if (initConditionParams.type != seissol::initializer::parameters::InitializationType::Zero &&
        !initConditionParams.avoidIC) {
      seissol::initializer::projectInitialField(
          memoryManager.initialConditions(), seissolInstance.meshReader(), storage);
    }
  }
}

void initSource(seissol::SeisSol& seissolInstance) {
  const auto& srcparams = seissolInstance.parameters().source;
  auto& memoryManager = seissolInstance.memoryManager();
  seissol::sourceterm::Manager::loadSources(srcparams.type,
                                            srcparams.fileName.c_str(),
                                            seissolInstance.meshReader(),
                                            memoryManager.ltsStorage(),
                                            memoryManager.backmap(),
                                            memoryManager.simulationSettings(),
                                            seissolInstance.timeManager());
}

} // namespace

void initSideConditions(seissol::SeisSol& seissolInstance) {
  logInfo() << "Setting initial conditions.";
  initInitialCondition(seissolInstance);
  logInfo() << "Reading source.";
  initSource(seissolInstance);
}

} // namespace seissol::initializer::initprocedure
