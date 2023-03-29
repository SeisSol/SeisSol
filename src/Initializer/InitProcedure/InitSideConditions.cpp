#include "SeisSol.h"

#include "Init.hpp"
#include "InitSideConditions.hpp"

#include "Initializer/InitialFieldProjection.h"
#include "Initializer/InputParameters.hpp"

#include "Parallel/MPI.h"

TravellingWaveParameters getTravellingWaveInformation() {
  const auto& icparams = seissol::SeisSol::main.getSeisSolParameters().initialization;

  TravellingWaveParameters m_travellingWaveParameters;
  m_travellingWaveParameters.origin = icparams.origin;
  m_travellingWaveParameters.kVec = icparams.kVec;
  constexpr double eps = 1e-15;
  for (size_t i = 0; i < NUMBER_OF_QUANTITIES; i++) {
    if (std::abs(m_travellingWaveParameters.ampField[i]) > eps) {
      m_travellingWaveParameters.varField.push_back(i);
      m_travellingWaveParameters.ampField.push_back(icparams.ampField[i]);
    }
  }
  return m_travellingWaveParameters;
}

std::vector<std::unique_ptr<physics::InitialField>> buildInitialConditionList() {
  const auto& icparams = seissol::SeisSol::main.getSeisSolParameters().initialization;
  auto& memmng = seissol::SeisSol::main.getMemoryManager();
  std::vector<std::unique_ptr<physics::InitialField>> iniconds;
  std::string initialConditionDescription = "";
  if (icparams.type == seissol::initializer::parameters::InitializationType::Planarwave) {
    initialConditionDescription = "Planar wave";
    auto materialData = memmng.getLtsLut()->lookup(memmng.getLts()->material, 0);

#ifdef MULTIPLE_SIMULATIONS
    for (int s = 0; s < MULTIPLE_SIMULATIONS; ++s) {
      const double phase = (2.0 * M_PI * s) / MULTIPLE_SIMULATIONS;
      iniconds.emplace_back(new physics::Planarwave(materialData, phase));
    }
#else
    iniconds.emplace_back(new physics::Planarwave(materialData));
#endif
  } else if (icparams.type ==
             seissol::initializer::parameters::InitializationType::SuperimposedPlanarwave) {
    initialConditionDescription = "Super-imposed planar wave";
#ifdef MULTIPLE_SIMULATIONS
    for (int s = 0; s < MULTIPLE_SIMULATIONS; ++s) {
      iniconds.emplace_back(new physics::SuperimposedPlanarwave(
          memmng.getLtsLut()->lookup(memmng.getLts()->material, 0),
          (2.0 * M_PI * s) / MULTIPLE_SIMULATIONS));
    }
#else
    iniconds.emplace_back(new physics::SuperimposedPlanarwave(
        memmng.getLtsLut()->lookup(memmng.getLts()->material, 0)));
#endif
  } else if (icparams.type == seissol::initializer::parameters::InitializationType::Zero) {
    initialConditionDescription = "Zero";
    iniconds.emplace_back(new physics::ZeroField());
  }
#if NUMBER_OF_RELAXATION_MECHANISMS == 0
  else if (icparams.type == seissol::initializer::parameters::InitializationType::Travelling) {
    initialConditionDescription = "Travelling wave";
    auto m_travellingWaveParameters = getTravellingWaveInformation();
    iniconds.emplace_back(new physics::TravellingWave(
        memmng.getLtsLut()->lookup(memmng.getLts()->material, 0), m_travellingWaveParameters));
  } else if (icparams.type == seissol::initializer::parameters::InitializationType::Scholte) {
    initialConditionDescription = "Scholte wave (elastic-acoustic)";
    iniconds.emplace_back(new physics::ScholteWave());
  } else if (icparams.type == seissol::initializer::parameters::InitializationType::Snell) {
    initialConditionDescription = "Snell's law (elastic-acoustic)";
    iniconds.emplace_back(new physics::SnellsLaw());
  } else if (icparams.type == seissol::initializer::parameters::InitializationType::Ocean0) {
    initialConditionDescription =
        "Ocean, an uncoupled ocean test case for acoustic equations (mode 0)";
    const auto g = seissol::SeisSol::main.getGravitationSetup().acceleration;
    iniconds.emplace_back(new physics::Ocean(0, g));
  } else if (icparams.type == seissol::initializer::parameters::InitializationType::Ocean1) {
    initialConditionDescription =
        "Ocean, an uncoupled ocean test case for acoustic equations (mode 1)";
    const auto g = seissol::SeisSol::main.getGravitationSetup().acceleration;
    iniconds.emplace_back(new physics::Ocean(1, g));
  } else if (icparams.type == seissol::initializer::parameters::InitializationType::Ocean2) {
    initialConditionDescription =
        "Ocean, an uncoupled ocean test case for acoustic equations (mode 2)";
    const auto g = seissol::SeisSol::main.getGravitationSetup().acceleration;
    iniconds.emplace_back(new physics::Ocean(2, g));
  }
#endif // NUMBER_OF_RELAXATION_MECHANISMS == 0
  else {
    logError() << "Non-implemented initial condition type:" << static_cast<int>(icparams.type);
  }
  logInfo(MPI::mpi.rank()) << "Using initial condition" << initialConditionDescription << ".";
  return iniconds;
}

void initInitialCondition() {
  auto iniconds = buildInitialConditionList();
  const auto& icparams = seissol::SeisSol::main.getSeisSolParameters().initialization;
  auto& memmng = seissol::SeisSol::main.getMemoryManager();

  if (icparams.type != seissol::initializer::parameters::InitializationType::Zero) {
    seissol::initializers::projectInitialField(iniconds,
                                               *memmng.getGlobalDataOnHost(),
                                               seissol::SeisSol::main.meshReader(),
                                               *memmng.getLts(),
                                               *memmng.getLtsLut());
  }

  memmng.setInitialConditions(std::move(iniconds));
}

void initSource() {
  const auto& srcparams = seissol::SeisSol::main.getSeisSolParameters().source;
  auto& memmng = seissol::SeisSol::main.getMemoryManager();
  SeisSol::main.sourceTermManager().loadSources(srcparams.type,
                                                srcparams.fileName.c_str(),
                                                seissol::SeisSol::main.meshReader(),
                                                memmng.getLtsTree(),
                                                memmng.getLts(),
                                                memmng.getLtsLut(),
                                                seissol::SeisSol::main.timeManager());
}

void initBoundary() {
  const auto& ssp = seissol::SeisSol::main.getSeisSolParameters();
  if (ssp.model.boundaryFileName !=
      "") { // TODO: better check than != ""; maybe some `enabled` field?
    seissol::SeisSol::main.getMemoryManager().initializeEasiBoundaryReader(
        ssp.model.boundaryFileName.c_str());
  }
}

void seissol::initializer::initprocedure::initSideConditions() {
  logInfo(seissol::MPI::mpi.rank()) << "Setting initial conditions.";
  initInitialCondition();
  logInfo(seissol::MPI::mpi.rank()) << "Reading source.";
  initSource();
  logInfo(seissol::MPI::mpi.rank()) << "Setting up boundary conditions.";
  initBoundary();
}
