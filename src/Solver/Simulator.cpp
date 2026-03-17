// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#include "Simulator.h"

#include "Memory/Tree/Layer.h"
#include "Modules/Modules.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Stopwatch.h"
#include "ResultWriter/AnalysisWriter.h"
#include "ResultWriter/EnergyOutput.h"
#include "SeisSol.h"
#include "TimeStepping/TimeManager.h"

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
  faultOutputManager->writePickpointOutput(0.0, 0.0);

  Stopwatch simulationStopwatch;
  if (checkpoint) {
    logInfo() << "The simulation will run from" << currentTime << "s (checkpoint time) until"
              << finalTime << "s.";
  } else {
    logInfo() << "The simulation will run until" << finalTime << "s.";
  }
  simulationStopwatch.start();

  Stopwatch computeStopwatch;
  Stopwatch ioStopwatch;

  ioStopwatch.start();

  // Set start time (required when loading a check)
  seissolInstance.timeManager().setInitialTimes(currentTime);

  const double timeTolerance = seissolInstance.timeManager().getTimeTolerance();

  // Write initial wave field snapshot
  if (checkpoint) {
    Modules::callSimulationStartHook(currentTime);
  } else {
    Modules::callSimulationStartHook(std::optional<double>{});
  }

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
      logError() << "Simulator did not advance in time from" << currentTime << "s to"
                 << upcomingTime << "s.";
    }
    if (aborted) {
      logInfo() << "Aborting simulation.";
      break;
    }

    // update the DOFs
    logInfo() << "Start simulation epoch. (from" << currentTime << "s to" << upcomingTime << "s)";
    computeStopwatch.start();
    seissolInstance.timeManager().advanceInTime(upcomingTime);
    computeStopwatch.pause();
    logInfo() << "End simulation epoch. (at" << upcomingTime << "s)";

    ioStopwatch.start();

    // update current time
    currentTime = upcomingTime;

    // Synchronize data (TODO(David): synchronize lazily)
    seissolInstance.timeManager().synchronizeTo(initializer::AllocationPlace::Host);

    // Check all synchronization point hooks and set the new upcoming time
    upcomingTime = std::min(finalTime, Modules::callSyncHook(currentTime, timeTolerance));

    ioStopwatch.pause();

    const double currentSplit = simulationStopwatch.split();
    Stopwatch::print("Time spent this epoch (total):", currentSplit - lastSplit);
    Stopwatch::print("Time spent this epoch (compute):", computeStopwatch.split());
    Stopwatch::print("Time spent this epoch (blocking IO):", ioStopwatch.split());
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

  const auto& outputParams = seissolInstance.getSeisSolParameters().output;

  const bool isLoopStatisticsNetcdfOutputOn = outputParams.loopStatisticsNetcdfOutput;
  const auto& outputPrefix = outputParams.prefix;
  seissolInstance.timeManager().printComputationTime(outputPrefix, isLoopStatisticsNetcdfOutputOn);

  seissolInstance.analysisWriter().printAnalysis(currentTime);

  seissolInstance.flopCounter().printPerformanceSummary(wallTime);
}

} // namespace seissol
