// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InitSideConditions.h"

#include "Equations/Datastructures.h"
#include "Initializer/InitialFieldProjection.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/LTS.h"
#include "Physics/InitialField.h"
#include "Physics/Scenario/Registry.h"
#include "SeisSol.h"
#include "SourceTerm/Manager.h"

#include <memory>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer::initprocedure {

namespace {

std::vector<std::unique_ptr<physics::InitialField>>
    buildInitialConditionList(seissol::SeisSol& seissolInstance) {
  const auto& parameters = seissolInstance.parameters();
  auto& memoryManager = seissolInstance.memoryManager();
  const auto type = parameters.initialization.type;

  const auto availability = physics::scenario::availability(type);
  if (!availability.available) {
    logError() << "The initial condition" << physics::scenario::name(type).data()
               << "cannot be used with material" << model::MaterialT::Text.c_str() << "--"
               << availability.reason.data() << ".";
  }

  const auto pos = memoryManager.backmap().get(0);
  const auto materialData = memoryManager.ltsStorage().lookup<LTS::Material>(pos);

  logInfo() << "Using initial condition" << physics::scenario::name(type).data() << ".";

  return physics::scenario::build(
      type, physics::scenario::Input{parameters, materialData, seissolInstance.gravitationSetup()});
}

void initInitialCondition(seissol::SeisSol& seissolInstance) {
  const auto& initConditionParams = seissolInstance.parameters().initialization;
  auto& memoryManager = seissolInstance.memoryManager();

  if (initConditionParams.type == seissol::initializer::parameters::InitializationType::Easi) {
    logInfo() << "Loading the initial condition from the easi file" << initConditionParams.filename;
    seissol::initializer::projectEasiInitialField({initConditionParams.filename},
                                                  *memoryManager.globalData().onHost,
                                                  seissolInstance.meshReader(),
                                                  memoryManager.ltsStorage(),
                                                  initConditionParams.hasTime);
  } else {
    auto initConditions = buildInitialConditionList(seissolInstance);
    if (initConditionParams.type != seissol::initializer::parameters::InitializationType::Zero &&
        !initConditionParams.avoidIC) {
      seissol::initializer::projectInitialField(initConditions,
                                                *memoryManager.globalData().onHost,
                                                seissolInstance.meshReader(),
                                                memoryManager.ltsStorage());
    }
    memoryManager.setInitialConditions(std::move(initConditions));
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
  logInfo() << "Reading source.";
  initSource(seissolInstance);
}

} // namespace seissol::initializer::initprocedure
