#include "SeisSol.h"

#include "Init.hpp"
#include "InitSideConditions.hpp"

#include "Initializer/InitialFieldProjection.h"
#include "Initializer/parameters/SeisSolParameters.h"

#include "Parallel/MPI.h"

namespace {

static TravellingWaveParameters getTravellingWaveInformation(seissol::SeisSol& seissolInstance) {
  const auto& initConditionParams = seissolInstance.getSeisSolParameters().initialization;

  TravellingWaveParameters travellingWaveParameters;
  travellingWaveParameters.origin = initConditionParams.origin;
  travellingWaveParameters.kVec = initConditionParams.kVec;
  constexpr double eps = 1e-15;
  for (size_t i = 0; i < NUMBER_OF_QUANTITIES; i++) {
    if (std::abs(initConditionParams.ampField[i]) > eps) {
      travellingWaveParameters.varField.push_back(i);
      travellingWaveParameters.ampField.push_back(initConditionParams.ampField[i]);
    }
  }
  return travellingWaveParameters;
}

static AcousticTravellingWaveParametersITM getAcousticTravellingWaveITMInformation(seissol::SeisSol& seissolInstance) {
  const auto& initConditionParams = seissolInstance.getSeisSolParameters().initialization;
  const auto& itmParams = seissolInstance.getSeisSolParameters().model.itmParameters;

  AcousticTravellingWaveParametersITM acousticTravellingWaveParametersITM;
  acousticTravellingWaveParametersITM.k = initConditionParams.k;
  acousticTravellingWaveParametersITM.itmStartingTime = itmParams.itmStartingTime;
  acousticTravellingWaveParametersITM.itmDuration = itmParams.itmDuration;
  acousticTravellingWaveParametersITM.itmVelocityScalingFactor = itmParams.itmVelocityScalingFactor;

  return acousticTravellingWaveParametersITM;
}

static std::vector<std::unique_ptr<physics::InitialField>>
    buildInitialConditionList(seissol::SeisSol& seissolInstance) {
  const auto& initConditionParams = seissolInstance.getSeisSolParameters().initialization;
  auto& memoryManager = seissolInstance.getMemoryManager();
  std::vector<std::unique_ptr<physics::InitialField>> initConditions;
  std::string initialConditionDescription = "";
  if (initConditionParams.type ==
      seissol::initializers::parameters::InitializationType::Planarwave) {
    initialConditionDescription = "Planar wave";
    auto materialData = memoryManager.getLtsLut()->lookup(memoryManager.getLts()->material, 0);

#ifdef MULTIPLE_SIMULATIONS
    for (int s = 0; s < MULTIPLE_SIMULATIONS; ++s) {
      const double phase = (2.0 * M_PI * s) / MULTIPLE_SIMULATIONS;
      initConditions.emplace_back(new physics::Planarwave(materialData, phase));
    }
#else
    initConditions.emplace_back(new physics::Planarwave(materialData));
#endif
  } else if (initConditionParams.type ==
             seissol::initializers::parameters::InitializationType::SuperimposedPlanarwave) {
    initialConditionDescription = "Super-imposed planar wave";
#ifdef MULTIPLE_SIMULATIONS
    for (int s = 0; s < MULTIPLE_SIMULATIONS; ++s) {
      initConditions.emplace_back(new physics::SuperimposedPlanarwave(
          memoryManager.getLtsLut()->lookup(memoryManager.getLts()->material, 0),
          (2.0 * M_PI * s) / MULTIPLE_SIMULATIONS));
    }
#else
    initConditions.emplace_back(new physics::SuperimposedPlanarwave(
        memoryManager.getLtsLut()->lookup(memoryManager.getLts()->material, 0)));
#endif
  } else if (initConditionParams.type ==
             seissol::initializers::parameters::InitializationType::Zero) {
    initialConditionDescription = "Zero";
    initConditions.emplace_back(new physics::ZeroField());
  }
#if NUMBER_OF_RELAXATION_MECHANISMS == 0
  else if (initConditionParams.type ==
           seissol::initializers::parameters::InitializationType::Travelling) {
    initialConditionDescription = "Travelling wave";
    auto travellingWaveParameters = getTravellingWaveInformation(seissolInstance);
    initConditions.emplace_back(new physics::TravellingWave(
        memoryManager.getLtsLut()->lookup(memoryManager.getLts()->material, 0),
        travellingWaveParameters));
  } else if (initConditionParams.type ==
             seissol::initializers::parameters::InitializationType::AcousticTravellingWithITM) {
    initialConditionDescription = "Acoustic Travelling Wave with ITM";
    auto acousticTravellingWaveParametersITM = getAcousticTravellingWaveITMInformation(seissolInstance);
    initConditions.emplace_back(new physics::AcousticTravellingWaveITM(
        memoryManager.getLtsLut()->lookup(memoryManager.getLts()->material, 0),
        acousticTravellingWaveParametersITM));
  } else if (initConditionParams.type ==
             seissol::initializers::parameters::InitializationType::Scholte) {
    initialConditionDescription = "Scholte wave (elastic-acoustic)";
    initConditions.emplace_back(new physics::ScholteWave());
  } else if (initConditionParams.type ==
             seissol::initializers::parameters::InitializationType::Snell) {
    initialConditionDescription = "Snell's law (elastic-acoustic)";
    initConditions.emplace_back(new physics::SnellsLaw());
  } else if (initConditionParams.type ==
             seissol::initializers::parameters::InitializationType::Ocean0) {
    initialConditionDescription =
        "Ocean, an uncoupled ocean test case for acoustic equations (mode 0)";
    const auto g = seissolInstance.getGravitationSetup().acceleration;
    initConditions.emplace_back(new physics::Ocean(0, g));
  } else if (initConditionParams.type ==
             seissol::initializers::parameters::InitializationType::Ocean1) {
    initialConditionDescription =
        "Ocean, an uncoupled ocean test case for acoustic equations (mode 1)";
    const auto g = seissolInstance.getGravitationSetup().acceleration;
    initConditions.emplace_back(new physics::Ocean(1, g));
  } else if (initConditionParams.type ==
             seissol::initializers::parameters::InitializationType::Ocean2) {
    initialConditionDescription =
        "Ocean, an uncoupled ocean test case for acoustic equations (mode 2)";
    const auto g = seissolInstance.getGravitationSetup().acceleration;
    initConditions.emplace_back(new physics::Ocean(2, g));
  } else if (initConditionParams.type ==
             seissol::initializers::parameters::InitializationType::PressureInjection) {
    initialConditionDescription = "Pressure Injection";
#ifndef USE_POROELASTIC
    logError()
        << "The initial condition 'Pressure Injection' only works with poroelastic materials.";
#endif
    initConditions.emplace_back(new physics::PressureInjection(initConditionParams));
  }
#endif // NUMBER_OF_RELAXATION_MECHANISMS == 0
  else {
    logError() << "Non-implemented initial condition type:"
               << static_cast<int>(initConditionParams.type);
  }
  logInfo(seissol::MPI::mpi.rank())
      << "Using initial condition" << initialConditionDescription << ".";
  return initConditions;
}

static void initInitialCondition(seissol::SeisSol& seissolInstance) {
  auto initConditions = buildInitialConditionList(seissolInstance);
  const auto& initConditionParams = seissolInstance.getSeisSolParameters().initialization;
  auto& memoryManager = seissolInstance.getMemoryManager();

  if (initConditionParams.type != seissol::initializers::parameters::InitializationType::Zero) {
    seissol::initializers::projectInitialField(initConditions,
                                               *memoryManager.getGlobalDataOnHost(),
                                               seissolInstance.meshReader(),
                                               seissolInstance.getMemoryManager(),
                                               *memoryManager.getLts(),
                                               *memoryManager.getLtsLut());
  }

  memoryManager.setInitialConditions(std::move(initConditions));
}

static void initSource(seissol::SeisSol& seissolInstance) {
  const auto& srcparams = seissolInstance.getSeisSolParameters().source;
  auto& memoryManager = seissolInstance.getMemoryManager();
  seissolInstance.sourceTermManager().loadSources(srcparams.type,
                                                  srcparams.fileName.c_str(),
                                                  seissolInstance.meshReader(),
                                                  memoryManager.getLtsTree(),
                                                  memoryManager.getLts(),
                                                  memoryManager.getLtsLut(),
                                                  seissolInstance.timeManager());
}

static void initBoundary(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  if (seissolParams.model.hasBoundaryFile) {
    seissolInstance.getMemoryManager().initializeEasiBoundaryReader(
        seissolParams.model.boundaryFileName.c_str());
  }
}

} // namespace

void seissol::initializers::initprocedure::initSideConditions(seissol::SeisSol& seissolInstance) {
  logInfo(seissol::MPI::mpi.rank()) << "Setting initial conditions.";
  initInitialCondition(seissolInstance);
  logInfo(seissol::MPI::mpi.rank()) << "Reading source.";
  initSource(seissolInstance);
  logInfo(seissol::MPI::mpi.rank()) << "Setting up boundary conditions.";
  initBoundary(seissolInstance);
}
