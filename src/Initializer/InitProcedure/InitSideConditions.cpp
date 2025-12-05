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

TravellingWaveParameters getTravellingWaveInformation(std::size_t quantities,
                                                      seissol::SeisSol& seissolInstance) {
  const auto& initConditionParams = seissolInstance.getSeisSolParameters().initialization;

  TravellingWaveParameters travellingWaveParameters{};
  travellingWaveParameters.origin = initConditionParams.origin;
  travellingWaveParameters.kVec = initConditionParams.kVec;
  constexpr double Eps = 1e-15;
  for (size_t i = 0; i < quantities; i++) {
    if (std::abs(initConditionParams.ampField[i]) > Eps) {
      travellingWaveParameters.varField.push_back(i);
      travellingWaveParameters.ampField.emplace_back(initConditionParams.ampField[i]);
    }
  }
  return travellingWaveParameters;
}

AcousticTravellingWaveParametersITM
    getAcousticTravellingWaveITMInformation(seissol::SeisSol& seissolInstance) {
  const auto& initConditionParams = seissolInstance.getSeisSolParameters().initialization;
  const auto& itmParams = seissolInstance.getSeisSolParameters().model.itmParameters;

  AcousticTravellingWaveParametersITM acousticTravellingWaveParametersITM{};
  acousticTravellingWaveParametersITM.k = initConditionParams.k;
  acousticTravellingWaveParametersITM.itmStartingTime = itmParams.itmStartingTime;
  acousticTravellingWaveParametersITM.itmDuration = itmParams.itmDuration;
  acousticTravellingWaveParametersITM.itmVelocityScalingFactor = itmParams.itmVelocityScalingFactor;

  return acousticTravellingWaveParametersITM;
}

std::vector<std::unique_ptr<physics::InitialField>>
    buildInitialConditionList(seissol::SeisSol& seissolInstance) {
  const auto& initConditionParams = seissolInstance.getSeisSolParameters().initialization;
  auto& memoryManager = seissolInstance.getMemoryManager();
  std::vector<std::unique_ptr<physics::InitialField>> initConditions;
  std::string initialConditionDescription;

  const auto pos = memoryManager.getBackmap().get(0);

  const auto& layer = memoryManager.getLtsStorage().layer(pos.color);

  layer.wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using MaterialT = model::MaterialTT<Cfg>;

    const auto& materialData = layer.var<LTS::Material>()[pos.cell];

    if (initConditionParams.type ==
        seissol::initializer::parameters::InitializationType::Planarwave) {
      initialConditionDescription = "Planar wave";

      for (std::size_t s = 0; s < seissol::multisim::NumSimulations<Cfg>; ++s) {
        const double phase = (2.0 * M_PI * s) / seissol::multisim::NumSimulations<Cfg>;
        initConditions.emplace_back(new physics::Planarwave<Cfg>(materialData, phase));
      }
    } else if (initConditionParams.type ==
               seissol::initializer::parameters::InitializationType::SuperimposedPlanarwave) {
      initialConditionDescription = "Super-imposed planar wave";

      for (std::size_t s = 0; s < seissol::multisim::NumSimulations<Cfg>; ++s) {
        const double phase = (2.0 * M_PI * s) / seissol::multisim::NumSimulations<Cfg>;
        initConditions.emplace_back(new physics::SuperimposedPlanarwave<Cfg>(materialData, phase));
      }
    } else if (initConditionParams.type ==
               seissol::initializer::parameters::InitializationType::Zero) {
      initialConditionDescription = "Zero";
      initConditions.emplace_back(new physics::ZeroField());
    } else if (initConditionParams.type ==
                   seissol::initializer::parameters::InitializationType::Travelling &&
               MaterialT::Mechanisms == 0) {
      initialConditionDescription = "Travelling wave";
      auto travellingWaveParameters =
          getTravellingWaveInformation(MaterialT::NumQuantities, seissolInstance);

      const auto materialData = memoryManager.getLtsStorage().lookup<LTS::Material>(pos);
      initConditions.emplace_back(
          new physics::TravellingWave<Cfg>(materialData, travellingWaveParameters));
    } else if (initConditionParams.type == seissol::initializer::parameters::InitializationType::
                                               AcousticTravellingWithITM &&
               MaterialT::Mechanisms == 0) {
      initialConditionDescription = "Acoustic Travelling Wave with ITM";
      auto acousticTravellingWaveParametersITM =
          getAcousticTravellingWaveITMInformation(seissolInstance);

      const auto materialData = memoryManager.getLtsStorage().lookup<LTS::Material>(pos);
      initConditions.emplace_back(new physics::AcousticTravellingWaveITM(
          materialData, acousticTravellingWaveParametersITM));
    } else if (initConditionParams.type ==
                   seissol::initializer::parameters::InitializationType::Scholte &&
               MaterialT::Mechanisms == 0) {
      initialConditionDescription = "Scholte wave (elastic-acoustic)";
      initConditions.emplace_back(new physics::ScholteWave());
    } else if (initConditionParams.type ==
                   seissol::initializer::parameters::InitializationType::Snell &&
               MaterialT::Mechanisms == 0) {
      initialConditionDescription = "Snell's law (elastic-acoustic)";
      initConditions.emplace_back(new physics::SnellsLaw());
    } else if (initConditionParams.type ==
                   seissol::initializer::parameters::InitializationType::Ocean0 &&
               MaterialT::Mechanisms == 0) {
      initialConditionDescription =
          "Ocean, an uncoupled ocean test case for acoustic equations (mode 0)";
      const auto g = seissolInstance.getGravitationSetup().acceleration;
      initConditions.emplace_back(new physics::Ocean(0, g));
    } else if (initConditionParams.type ==
                   seissol::initializer::parameters::InitializationType::Ocean1 &&
               MaterialT::Mechanisms == 0) {
      initialConditionDescription =
          "Ocean, an uncoupled ocean test case for acoustic equations (mode 1)";
      const auto g = seissolInstance.getGravitationSetup().acceleration;
      initConditions.emplace_back(new physics::Ocean(1, g));
    } else if (initConditionParams.type ==
                   seissol::initializer::parameters::InitializationType::Ocean2 &&
               MaterialT::Mechanisms == 0) {
      initialConditionDescription =
          "Ocean, an uncoupled ocean test case for acoustic equations (mode 2)";
      const auto g = seissolInstance.getGravitationSetup().acceleration;
      initConditions.emplace_back(new physics::Ocean(2, g));
    } else if (initConditionParams.type ==
                   seissol::initializer::parameters::InitializationType::PressureInjection &&
               MaterialT::Type == seissol::model::MaterialType::Poroelastic) {
      initialConditionDescription = "Pressure Injection";
      initConditions.emplace_back(new physics::PressureInjection(initConditionParams));
    } else {
      logError() << "Non-implemented initial condition type:"
                 << static_cast<int>(initConditionParams.type);
    }
    logInfo() << "Using initial condition" << initialConditionDescription << ".";
  });
  return initConditions;
}

void initInitialCondition(seissol::SeisSol& seissolInstance) {
  const auto& initConditionParams = seissolInstance.getSeisSolParameters().initialization;
  auto& memoryManager = seissolInstance.getMemoryManager();

  if (initConditionParams.type == seissol::initializer::parameters::InitializationType::Easi) {
    logInfo() << "Loading the initial condition from the easi file" << initConditionParams.filename;
    seissol::initializer::projectEasiInitialField({initConditionParams.filename},
                                                  memoryManager.getGlobalData(),
                                                  seissolInstance.meshReader(),
                                                  memoryManager.getLtsStorage(),
                                                  initConditionParams.hasTime);
  } else {
    auto initConditions = buildInitialConditionList(seissolInstance);
    if (initConditionParams.type != seissol::initializer::parameters::InitializationType::Zero &&
        !initConditionParams.avoidIC) {
      seissol::initializer::projectInitialField(initConditions,
                                                memoryManager.getGlobalData(),
                                                seissolInstance.meshReader(),
                                                memoryManager.getLtsStorage());
    }
    memoryManager.setInitialConditions(std::move(initConditions));
  }
}

void initSource(seissol::SeisSol& seissolInstance) {
  const auto& srcparams = seissolInstance.getSeisSolParameters().source;
  auto& memoryManager = seissolInstance.getMemoryManager();
  seissol::sourceterm::Manager::loadSources(srcparams.type,
                                            srcparams.fileName.c_str(),
                                            seissolInstance.meshReader(),
                                            memoryManager.getLtsStorage(),
                                            memoryManager.getBackmap(),
                                            seissolInstance.timeManager());
}

void initBoundary(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  if (seissolParams.model.hasBoundaryFile) {
    seissolInstance.getMemoryManager().initializeEasiBoundaryReader(
        seissolParams.model.boundaryFileName.c_str());
  }
}

} // namespace

void initSideConditions(seissol::SeisSol& seissolInstance) {
  logInfo() << "Setting initial conditions.";
  initInitialCondition(seissolInstance);
  logInfo() << "Reading source.";
  initSource(seissolInstance);
  logInfo() << "Setting up boundary conditions.";
  initBoundary(seissolInstance);
}

} // namespace seissol::initializer::initprocedure
