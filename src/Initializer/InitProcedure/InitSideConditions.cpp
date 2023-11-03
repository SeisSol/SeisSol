#include "SeisSol.h"

#include "Init.hpp"
#include "InitSideConditions.hpp"

#include "Initializer/InitialFieldProjection.h"
#include "Initializer/InputParameters.hpp"

#include "Parallel/MPI.h"

namespace {

static TravellingWaveParameters getTravellingWaveInformation() {
  const auto& initConditionParams = seissol::SeisSol::main.getSeisSolParameters().initialization;

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

static AcousticTravellingWaveParametersITM getAcousticTravellingWaveITMInformation(){
  const auto& initConditionParams = seissol::SeisSol::main.getSeisSolParameters().initialization;
        AcousticTravellingWaveParametersITM acousticTravellingWaveParametersITM;
        acousticTravellingWaveParametersITM.k = initConditionParams.k;

        return acousticTravellingWaveParametersITM;
}

static std::vector<std::unique_ptr<physics::InitialField>> buildInitialConditionList() {
  const auto& initConditionParams = seissol::SeisSol::main.getSeisSolParameters().initialization;
  auto& memoryManager = seissol::SeisSol::main.getMemoryManager();
  std::vector<std::unique_ptr<physics::InitialField>> initConditions;
  std::string initialConditionDescription = "";
  if (initConditionParams.type ==
      seissol::initializer::parameters::InitializationType::Planarwave) {
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
             seissol::initializer::parameters::InitializationType::SuperimposedPlanarwave) {
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
             seissol::initializer::parameters::InitializationType::Zero) {
    initialConditionDescription = "Zero";
    initConditions.emplace_back(new physics::ZeroField());
  }
#if NUMBER_OF_RELAXATION_MECHANISMS == 0
  else if (initConditionParams.type ==
           seissol::initializer::parameters::InitializationType::Travelling) {
    initialConditionDescription = "Travelling wave";
    auto travellingWaveParameters = getTravellingWaveInformation();
    initConditions.emplace_back(new physics::TravellingWave(
        memoryManager.getLtsLut()->lookup(memoryManager.getLts()->material, 0),
        travellingWaveParameters));
  } else if (initConditionParams.type == seissol::initializer::parameters::InitializationType::AcousticTravellingwithITM){
        initialConditionDescription = "Acoustic Travelling Wave with ITM";
        auto acousticTravellingWaveParametersITM = getAcousticTravellingWaveITMInformation();
        initConditions.emplace_back(new physics::AcousticTravellingWaveITM(
                memoryManager.getLtsLut()->lookup(memoryManager.getLts()->material, 0),
                acousticTravellingWaveParametersITM));
  } else if (initConditionParams.type ==
             seissol::initializer::parameters::InitializationType::Scholte) {
    initialConditionDescription = "Scholte wave (elastic-acoustic)";
    initConditions.emplace_back(new physics::ScholteWave());
  } else if (initConditionParams.type ==
             seissol::initializer::parameters::InitializationType::Snell) {
    initialConditionDescription = "Snell's law (elastic-acoustic)";
    initConditions.emplace_back(new physics::SnellsLaw());
  } else if (initConditionParams.type ==
             seissol::initializer::parameters::InitializationType::Ocean0) {
    initialConditionDescription =
        "Ocean, an uncoupled ocean test case for acoustic equations (mode 0)";
    const auto g = seissol::SeisSol::main.getGravitationSetup().acceleration;
    initConditions.emplace_back(new physics::Ocean(0, g));
  } else if (initConditionParams.type ==
             seissol::initializer::parameters::InitializationType::Ocean1) {
    initialConditionDescription =
        "Ocean, an uncoupled ocean test case for acoustic equations (mode 1)";
    const auto g = seissol::SeisSol::main.getGravitationSetup().acceleration;
    initConditions.emplace_back(new physics::Ocean(1, g));
  } else if (initConditionParams.type ==
             seissol::initializer::parameters::InitializationType::Ocean2) {
    initialConditionDescription =
        "Ocean, an uncoupled ocean test case for acoustic equations (mode 2)";
    const auto g = seissol::SeisSol::main.getGravitationSetup().acceleration;
    initConditions.emplace_back(new physics::Ocean(2, g));
  } else if (initConditionParams.type ==
             seissol::initializer::parameters::InitializationType::PressureInjection) {
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

static void initInitialCondition() {
  auto initConditions = buildInitialConditionList();
  const auto& initConditionParams = seissol::SeisSol::main.getSeisSolParameters().initialization;
  auto& memoryManager = seissol::SeisSol::main.getMemoryManager();

  if (initConditionParams.type != seissol::initializer::parameters::InitializationType::Zero) {
    seissol::initializers::projectInitialField(initConditions,
                                               *memoryManager.getGlobalDataOnHost(),
                                               seissol::SeisSol::main.meshReader(),
                                               *memoryManager.getLts(),
                                               *memoryManager.getLtsLut());
  }

  memoryManager.setInitialConditions(std::move(initConditions));
}

static void initSource() {
  const auto& srcparams = seissol::SeisSol::main.getSeisSolParameters().source;
  auto& memoryManager = seissol::SeisSol::main.getMemoryManager();
  SeisSol::main.sourceTermManager().loadSources(srcparams.type,
                                                srcparams.fileName.c_str(),
                                                seissol::SeisSol::main.meshReader(),
                                                memoryManager.getLtsTree(),
                                                memoryManager.getLts(),
                                                memoryManager.getLtsLut(),
                                                seissol::SeisSol::main.timeManager());
}

static void initBoundary() {
  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();
  if (seissolParams.model.hasBoundaryFile) {
    seissol::SeisSol::main.getMemoryManager().initializeEasiBoundaryReader(
        seissolParams.model.boundaryFileName.c_str());
  }
}

} // namespace

void seissol::initializer::initprocedure::initSideConditions() {
  logInfo(seissol::MPI::mpi.rank()) << "Setting initial conditions.";
  initInitialCondition();
  logInfo(seissol::MPI::mpi.rank()) << "Reading source.";
  initSource();
  logInfo(seissol::MPI::mpi.rank()) << "Setting up boundary conditions.";
  initBoundary();
}
