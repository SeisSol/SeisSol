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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer::initprocedure {

namespace {

/// Whether the scenario takes its medium from the material of a cell, as
/// opposed to stating its own or not needing one.
bool builtFromCellMaterial(parameters::InitializationType type) {
  switch (type) {
  case parameters::InitializationType::Planarwave:
  case parameters::InitializationType::SuperimposedPlanarwave:
  case parameters::InitializationType::Travelling:
  case parameters::InitializationType::AcousticTravellingWithITM:
    return true;
  default:
    return false;
  }
}

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

  // The scenarios built from the material of a cell are solutions of a
  // homogeneous medium, and they take that medium from one cell. Where the
  // material is allowed to vary inside a cell, that assumption is worth
  // checking: the samples of this cell have to agree with each other, or the
  // field is a solution of a medium that is not the one being simulated. This
  // says nothing about the other cells, which the choice of a single cell has
  // always assumed.
  if constexpr (NodalMaterial) {
    if (builtFromCellMaterial(type)) {
      const auto& samples = memoryManager.ltsStorage().lookup<LTS::NodalMaterialData>(pos);
      for (const auto& [name, member] : model::MaterialT::ParameterMap) {
        const double reference = samples[0].*member;
        const double scale = std::max(1.0, std::abs(reference));
        for (std::size_t node = 1; node < LTS::MaterialNodes; ++node) {
          if (std::abs(samples[node].*member - reference) > 1.0e-12 * scale) {
            logError() << "The initial condition" << physics::scenario::name(type).data()
                       << "is a solution of a homogeneous medium, but the material"
                       << "varies inside the cell it was built from (" << name.c_str()
                       << "). Use a material without sub-cell variation for this initial"
                       << "condition.";
          }
        }
      }
    }
  }

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
