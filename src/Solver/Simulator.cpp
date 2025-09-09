// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#include "Simulator.h"

#include "Modules/Modules.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Stopwatch.h"
#include "ResultWriter/AnalysisWriter.h"
#include "ResultWriter/EnergyOutput.h"
#include "SeisSol.h"
#include "TimeStepping/TimeManager.h"
#include <Memory/Tree/Layer.h>
#include <Monitoring/Instrumentation.h>
#include <Parallel/Runtime/Stream.h>
#include <algorithm>
#include <cassert>
#include <optional>
#include <utils/logger.h>

namespace seissol {

Simulator::Simulator() = default;

void Simulator::setFinalTime(double finalTime) {
  assert(finalTime > 0);
  this->finalTime = finalTime;
}

void Simulator::setUsePlasticity(bool plasticity) { usePlasticity = plasticity; }

void Simulator::setCurrentTime(double currentTime) {
  assert(currentTime >= 0);
  this->currentTime = currentTime;
  checkpoint = true;
}

void Simulator::abort() { aborted = true; }

void Simulator::simulate(SeisSol& seissolInstance) {
  SCOREP_USER_REGION("simulate", SCOREP_USER_REGION_TYPE_FUNCTION)

  auto* faultOutputManager = seissolInstance.timeManager().getFaultOutputManager();
  parallel::runtime::StreamRuntime runtime;
  faultOutputManager->writePickpointOutput(0.0, 0.0, runtime);
  runtime.wait();

  Stopwatch simulationStopwatch;
  simulationStopwatch.start();

  Stopwatch computeStopwatch;
  Stopwatch ioStopwatch;

  ioStopwatch.start();

  // Set start time (required for checkpointing)
  seissolInstance.timeManager().setInitialTimes(currentTime);

  const double timeTolerance = seissolInstance.timeManager().getTimeTolerance();

  // Write initial wave field snapshot
  if (checkpoint) {
    Modules::callSimulationStartHook(currentTime);
  } else {
    Modules::callSimulationStartHook(std::optional<double>{});
  }

  // intialize wave field and checkpoint time
  Modules::setSimulationStartTime(currentTime);

  // derive next synchronization time
  double upcomingTime = finalTime;
  // NOTE: This will not call the module specific implementation of the synchronization hook
  // since the current time is the simulation start time. We only use this function here to
  // get correct upcoming time. To be on the safe side, we use zero time tolerance.
  upcomingTime = std::min(upcomingTime, Modules::callSyncHook(currentTime, 0.0));

  double lastSplit = 0;

  ioStopwatch.pause();

  Stopwatch::print("Time spent for initial IO:", ioStopwatch.split());

  // use an empty log message as visual separator
  logInfo() << "";

  while (finalTime > currentTime + timeTolerance) {
    if (upcomingTime < currentTime + timeTolerance) {
      logError() << "Simulator did not advance in time from" << currentTime << "to" << upcomingTime;
    }
    if (aborted) {
      logInfo() << "Aborting simulation.";
      break;
    }

    // update the DOFs
    logInfo() << "Start simulation epoch.";
    computeStopwatch.start();
    seissolInstance.timeManager().advanceInTime(upcomingTime);
    computeStopwatch.pause();
    logInfo() << "End simulation epoch. Sync point.";

    ioStopwatch.start();

    // update current time
    currentTime = upcomingTime;

    // Synchronize data (TODO(David): synchronize lazily)
    seissolInstance.timeManager().synchronizeTo(initializer::AllocationPlace::Host);

    // Check all synchronization point hooks and set the new upcoming time
    upcomingTime = std::min(finalTime, Modules::callSyncHook(currentTime, timeTolerance));

    ioStopwatch.pause();

    const double currentSplit = simulationStopwatch.split();
    Stopwatch::print("Time spent this phase (total):", currentSplit - lastSplit);
    Stopwatch::print("Time spent this phase (compute):", computeStopwatch.split());
    Stopwatch::print("Time spent this phase (blocking IO):", ioStopwatch.split());
    seissolInstance.flopCounter().printPerformanceUpdate(currentSplit);
    lastSplit = currentSplit;

    // use an empty log message as visual separator
    logInfo() << "";
  }

  // synchronize data (TODO(David): synchronize lazily)
  seissolInstance.timeManager().synchronizeTo(initializer::AllocationPlace::Host);

  Modules::callSyncHook(currentTime, timeTolerance, true);

  const double wallTime = simulationStopwatch.pause();
  simulationStopwatch.printTime("Simulation time (total):");
  computeStopwatch.printTime("Simulation time (compute):");
  ioStopwatch.printTime("Simulation time (blocking IO):");

  Modules::callHook<ModuleHook::SimulationEnd>();

  const auto& memoryManager = seissolInstance.getMemoryManager();
  const bool isLoopStatisticsNetcdfOutputOn = memoryManager.isLoopStatisticsNetcdfOutputOn();
  const auto& outputPrefix = memoryManager.getOutputPrefix();
  seissolInstance.timeManager().printComputationTime(outputPrefix, isLoopStatisticsNetcdfOutputOn);

  seissolInstance.analysisWriter().printAnalysis(currentTime);

  seissolInstance.flopCounter().printPerformanceSummary(wallTime);
}

} // namespace seissol
