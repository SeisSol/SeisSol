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
#include "Model/CommonDatastructures.h"
#include "Physics/InitialField.h"
#include "SeisSol.h"
#include "Solver/MultipleSimulations.h"
#include "SourceTerm/Manager.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <math.h>
#include <memory>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer::initprocedure {

namespace {

TravellingWaveParameters getTravellingWaveInformation(seissol::SeisSol& seissolInstance) {
  const auto& initConditionParams = seissolInstance.parameters().initialization;

  TravellingWaveParameters travellingWaveParameters{};
  travellingWaveParameters.origin = initConditionParams.origin;
  travellingWaveParameters.kVec = initConditionParams.kVec;
  constexpr double Eps = 1e-15;
  for (size_t i = 0; i < seissol::model::MaterialT::NumQuantities; i++) {
    if (std::abs(initConditionParams.ampField[i]) > Eps) {
      travellingWaveParameters.varField.push_back(i);
      travellingWaveParameters.ampField.emplace_back(initConditionParams.ampField[i]);
    }
  }
  return travellingWaveParameters;
}

AcousticTravellingWaveParametersITM
    getAcousticTravellingWaveITMInformation(seissol::SeisSol& seissolInstance) {
  const auto& initConditionParams = seissolInstance.parameters().initialization;
  const auto& itmParams = seissolInstance.parameters().model.itmParameters;

  AcousticTravellingWaveParametersITM acousticTravellingWaveParametersITM{};
  acousticTravellingWaveParametersITM.k = initConditionParams.k;
  acousticTravellingWaveParametersITM.itmStartingTime = itmParams.itmStartingTime;
  acousticTravellingWaveParametersITM.itmDuration = itmParams.itmDuration;
  acousticTravellingWaveParametersITM.itmVelocityScalingFactor = itmParams.itmVelocityScalingFactor;

  return acousticTravellingWaveParametersITM;
}

std::vector<std::unique_ptr<physics::InitialField>>
    buildInitialConditionList(seissol::SeisSol& seissolInstance) {
  const auto& initConditionParams = seissolInstance.parameters().initialization;
  auto& memoryManager = seissolInstance.memoryManager();
  std::vector<std::unique_ptr<physics::InitialField>> initConditions;
  std::string initialConditionDescription;

  const auto pos = memoryManager.backmap().get(0);

  // The analytical fields below are solutions of a homogeneous medium, and every
  // one of them is built from the material of one cell. Where the material is
  // allowed to vary inside a cell, that assumption is worth checking: the
  // samples of this cell have to agree with each other, or the field is a
  // solution of a medium that is not the one being simulated. This says nothing
  // about the other cells, which the choice of a single cell has always assumed.
  const auto homogeneousOrFail = [&](const std::string& description) {
    if constexpr (NodalMaterial) {
      const auto& samples = memoryManager.ltsStorage().lookup<LTS::NodalMaterialData>(pos);
      for (const auto& [name, member] : model::MaterialT::ParameterMap) {
        const double reference = samples[0].*member;
        const double scale = std::max(1.0, std::abs(reference));
        for (std::size_t node = 1; node < LTS::MaterialNodes; ++node) {
          if (std::abs(samples[node].*member - reference) > 1.0e-12 * scale) {
            logError() << description << "is a solution of a homogeneous medium, but the material"
                       << "varies inside the cell it was built from (" << name.c_str()
                       << "). Use a material without sub-cell variation for this initial"
                       << "condition.";
          }
        }
      }
    }
  };

  if (initConditionParams.type ==
      seissol::initializer::parameters::InitializationType::Planarwave) {
    initialConditionDescription = "Planar wave";
    homogeneousOrFail(initialConditionDescription);
    const auto materialData = memoryManager.ltsStorage().lookup<LTS::Material>(pos);

    for (std::size_t s = 0; s < seissol::multisim::NumSimulations; ++s) {
      const double phase = (2.0 * M_PI * s) / seissol::multisim::NumSimulations;
      initConditions.emplace_back(new physics::Planarwave(materialData, phase));
    }
  } else if (initConditionParams.type ==
             seissol::initializer::parameters::InitializationType::SuperimposedPlanarwave) {
    initialConditionDescription = "Super-imposed planar wave";
    homogeneousOrFail(initialConditionDescription);

    const auto materialData = memoryManager.ltsStorage().lookup<LTS::Material>(pos);
    for (std::size_t s = 0; s < seissol::multisim::NumSimulations; ++s) {
      const double phase = (2.0 * M_PI * s) / seissol::multisim::NumSimulations;
      initConditions.emplace_back(new physics::SuperimposedPlanarwave(materialData, phase));
    }
  } else if (initConditionParams.type ==
             seissol::initializer::parameters::InitializationType::Zero) {
    initialConditionDescription = "Zero";
    initConditions.emplace_back(new physics::ZeroField());
  } else if (initConditionParams.type ==
                 seissol::initializer::parameters::InitializationType::Travelling &&
             model::MaterialT::Mechanisms == 0) {
    initialConditionDescription = "Travelling wave";
    homogeneousOrFail(initialConditionDescription);
    auto travellingWaveParameters = getTravellingWaveInformation(seissolInstance);

    const auto materialData = memoryManager.ltsStorage().lookup<LTS::Material>(pos);
    initConditions.emplace_back(
        new physics::TravellingWave(materialData, travellingWaveParameters));
  } else if (initConditionParams.type ==
                 seissol::initializer::parameters::InitializationType::AcousticTravellingWithITM &&
             model::MaterialT::Mechanisms == 0) {
    initialConditionDescription = "Acoustic Travelling Wave with ITM";
    homogeneousOrFail(initialConditionDescription);
    auto acousticTravellingWaveParametersITM =
        getAcousticTravellingWaveITMInformation(seissolInstance);

    const auto materialData = memoryManager.ltsStorage().lookup<LTS::Material>(pos);
    initConditions.emplace_back(
        new physics::AcousticTravellingWaveITM(materialData, acousticTravellingWaveParametersITM));
  } else if (initConditionParams.type ==
                 seissol::initializer::parameters::InitializationType::Scholte &&
             model::MaterialT::Mechanisms == 0) {
    initialConditionDescription = "Scholte wave (elastic-acoustic)";
    initConditions.emplace_back(new physics::ScholteWave());
  } else if (initConditionParams.type ==
                 seissol::initializer::parameters::InitializationType::Snell &&
             model::MaterialT::Mechanisms == 0) {
    initialConditionDescription = "Snell's law (elastic-acoustic)";
    initConditions.emplace_back(new physics::SnellsLaw());
  } else if (initConditionParams.type ==
                 seissol::initializer::parameters::InitializationType::Ocean0 &&
             model::MaterialT::Mechanisms == 0) {
    initialConditionDescription =
        "Ocean, an uncoupled ocean test case for acoustic equations (mode 0)";
    const auto g = seissolInstance.gravitationSetup().acceleration;
    initConditions.emplace_back(new physics::Ocean(0, g));
  } else if (initConditionParams.type ==
                 seissol::initializer::parameters::InitializationType::Ocean1 &&
             model::MaterialT::Mechanisms == 0) {
    initialConditionDescription =
        "Ocean, an uncoupled ocean test case for acoustic equations (mode 1)";
    const auto g = seissolInstance.gravitationSetup().acceleration;
    initConditions.emplace_back(new physics::Ocean(1, g));
  } else if (initConditionParams.type ==
                 seissol::initializer::parameters::InitializationType::Ocean2 &&
             model::MaterialT::Mechanisms == 0) {
    initialConditionDescription =
        "Ocean, an uncoupled ocean test case for acoustic equations (mode 2)";
    const auto g = seissolInstance.gravitationSetup().acceleration;
    initConditions.emplace_back(new physics::Ocean(2, g));
  } else if (initConditionParams.type ==
                 seissol::initializer::parameters::InitializationType::PressureInjection &&
             model::MaterialT::Type == model::MaterialType::Poroelastic) {
    initialConditionDescription = "Pressure Injection";
    initConditions.emplace_back(new physics::PressureInjection(initConditionParams));
  } else {
    logError() << "Non-implemented initial condition type:"
               << static_cast<int>(initConditionParams.type);
  }
  logInfo() << "Using initial condition" << initialConditionDescription << ".";
  return initConditions;
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
