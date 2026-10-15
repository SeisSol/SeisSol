// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitSideConditions.h"

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Equations/Datastructures.h"
#include "Expr/RtcGpuDriver.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/InitialFieldProjection.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/ScriptParameters.h"
#include "Model/StateNames.h"
#include "Physics/GhostState.h"
#include "Physics/InitialField.h"
#include "Physics/Scenario/Registry.h"
#include "Physics/ScriptField.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"
#include "SourceTerm/Manager.h"

#include <algorithm>
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
/// `materialData` of one of them.
std::vector<std::unique_ptr<physics::InitialField>> buildInitialConditionList(
    seissol::SeisSol& seissolInstance, ConfigId config, const CellMaterialData& materialData) {
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
      physics::scenario::Input{
          parameters, materialData, seissolInstance.gravitationSetup(), config});
}

void initInitialCondition(seissol::SeisSol& seissolInstance) {
  const auto& initConditionParams = seissolInstance.parameters().initialization;
  auto& memoryManager = seissolInstance.memoryManager();

  if (initConditionParams.type == seissol::initializer::parameters::InitializationType::Script) {
    logInfo() << "Loading the initial condition from the script" << initConditionParams.filename;
    seissol::initializer::projectScriptInitialField(initConditionParams.filename,
                                                    seissolInstance.meshReader(),
                                                    memoryManager.ltsStorage(),
                                                    initConditionParams.hasTime);

    // The same script, as the field an analytic boundary condition asks for at its points and
    // times, one per fused simulation.
    for (const auto config : seissolInstance.parameters().model.configs()) {
      dispatchConfig(config, [&](auto cfg) {
        using Cfg = decltype(cfg);
        std::vector<std::unique_ptr<physics::InitialField>> fields;
        for (std::size_t sim = 0; sim < Cfg::NumSimulations; ++sim) {
          fields.push_back(std::make_unique<physics::ScriptField>(
              initConditionParams.filename,
              model::stateNames<model::MaterialOf<Cfg>>(),
              model::scriptParameterNames<model::MaterialOf<Cfg>>(),
              sim,
              initConditionParams.hasTime));
        }
        memoryManager.setInitialConditions(config, std::move(fields));
      });
    }
  } else {
    logInfo() << "Using initial condition"
              << physics::scenario::name(initConditionParams.type).data() << ".";

    // The configurations of one material set up the scenario alike, with the material of the
    // first cell of the first of them in the run that has cells on this rank: setting up a
    // scenario may pick among degenerate eigenvectors of the material, which then roundoff
    // decides, and configurations can average the material of a cell differently (e.g. with the
    // quadrature of another order).
    auto& storage = memoryManager.ltsStorage();
    auto& backmap = memoryManager.backmap();
    const auto cellCount = seissolInstance.meshReader().getElements().size();
    std::vector<std::optional<std::size_t>> firstCell(builtConfigCount());
    for (std::size_t cell = 0; cell < cellCount; ++cell) {
      const auto config = storage.lookup<LTS::SecondaryInformation>(backmap.get(cell)).configId;
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
              buildInitialConditionList(
                  seissolInstance,
                  config,
                  storage.lookup<LTS::Material>(backmap.get(firstCell[reference].value()))));
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

/// The ghost state of the analytical boundary integrated over a step at once, for every
/// configuration with such a face whose initial condition is a script that compiles.
void initAnalyticalBoundary(seissol::SeisSol& seissolInstance) {
  auto& memoryManager = seissolInstance.memoryManager();

  std::vector<bool> present(builtConfigCount(), false);
  for (auto& layer : memoryManager.ltsStorage().leaves(Ghost)) {
    const auto* cellInformation = layer.var<LTS::CellInformation>();
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      for (const auto faceType : cellInformation[cell].faceTypes) {
        if (faceType == FaceType::Analytical) {
          present[layer.getIdentifier().config] = true;
        }
      }
    }
  }

  for (const auto config : seissolInstance.parameters().model.configs()) {
    if (!present[config]) {
      continue;
    }
    dispatchConfig(config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      // the components of the state at a face, which leaves out the memory variables that the
      // solver keeps apart
      auto components = model::stateNames<model::MaterialOf<Cfg>>();
      components.resize(std::min<std::size_t>(
          components.size(), tensor::INodal<Cfg>::Shape[multisim::BasisDim<Cfg> + 1]));
      auto state = physics::AnalyticalGhostState::fromFields(
          memoryManager.initialConditions(config),
          components,
          model::scriptParameterNames<model::MaterialOf<Cfg>>(),
          Cfg::NumSimulations,
          Cfg::ConvergenceOrder);
      if constexpr (isDeviceOn()) {
        // a device build evaluates the boundary on the device, and nowhere else
        if (state == nullptr) {
          logError() << "The analytical boundary of a device build evaluates the initial "
                        "condition on the device, as a script that compiles: an sderiv module, a "
                        "Lua model that traces, or a scenario of SeisSol. An easi file or a Lua "
                        "model that does not trace is none.";
        }
        if (!expr::deviceRuntimeAvailable()) {
          logError() << "The analytical boundary of a device build needs the runtime compiler of "
                        "the device, which this build does not have.";
        }
      }
      memoryManager.setAnalyticalGhostState(config, std::move(state));
    });
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
                                            seissolInstance.timeManager());
}

} // namespace

void initSideConditions(seissol::SeisSol& seissolInstance) {
  logInfo() << "Setting initial conditions.";
  initInitialCondition(seissolInstance);
  initAnalyticalBoundary(seissolInstance);
  logInfo() << "Reading source.";
  initSource(seissolInstance);
}

} // namespace seissol::initializer::initprocedure
