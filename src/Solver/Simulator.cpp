// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 **/

#include <Initializer/Tree/Layer.h>
#include <Parallel/Runtime/Stream.h>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <limits>
#include <utils/logger.h>
#include <xdmfwriter/scorep_wrapper.h>

#include "Clustering/TimeManager.h"
#include "Modules/Modules.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Stopwatch.h"
#include "Monitoring/Unit.h"
#include "ResultWriter/AnalysisWriter.h"
#include "ResultWriter/EnergyOutput.h"
#include "SeisSol.h"
#include "Simulator.h"

seissol::Simulator::Simulator()
    : currentTime(0), finalTime(0), usePlasticity(false), checkPointTime(0),
      checkPointInterval(std::numeric_limits<double>::max()), aborted(false) {}

void seissol::Simulator::setFinalTime(double finalTime) {
  assert(finalTime > 0);
  this->finalTime = finalTime;
}

void seissol::Simulator::setUsePlasticity(bool plasticity) { usePlasticity = plasticity; }

void seissol::Simulator::setCurrentTime(double currentTime) {
  assert(currentTime > 0);
  this->currentTime = currentTime;
}

void seissol::Simulator::abort() { aborted = true; }

void seissol::Simulator::simulate(seissol::SeisSol& seissolInstance) {
  SCOREP_USER_REGION("simulate", SCOREP_USER_REGION_TYPE_FUNCTION)

  auto* faultOutputManager = seissolInstance.timeManager().getFaultOutputManager();
  parallel::runtime::StreamRuntime syncRuntime;
  faultOutputManager->writePickpointOutput(0.0, 0.0, syncRuntime);
  syncRuntime.wait();

  Stopwatch simulationStopwatch;
  simulationStopwatch.start();

  Stopwatch computeStopwatch;
  Stopwatch ioStopwatch;

  ioStopwatch.start();

  // Set start time (required for checkpointing)
  seissolInstance.timeManager().setInitialTimes(currentTime);

  const double timeTolerance = seissolInstance.timeManager().getTimeTolerance();

  // Write initial wave field snapshot
  if (currentTime == 0.0) {
    Modules::callHook<ModuleHook::SimulationStart>();
  }

  // intialize wave field and checkpoint time
  Modules::setSimulationStartTime(currentTime);

  // derive next synchronization time
  double upcomingTime = finalTime;
  // NOTE: This will not call the module specific implementation of the synchronization hook
  // since the current time is the simulation start time. We only use this function here to
  // get correct upcoming time. To be on the safe side, we use zero time tolerance.
  upcomingTime = std::min( upcomingTime, Modules::callSyncHook(currentTime, 0.0) );

  double lastSplit = 0;

  ioStopwatch.pause();

  Stopwatch::print("Time spent for initial IO:", ioStopwatch.split(), seissol::MPI::mpi.comm());

  while (finalTime > currentTime + timeTolerance) {
    if (upcomingTime < currentTime + timeTolerance) {
      logError() << "Simulator did not advance in time from" << currentTime << "to" << upcomingTime;
    }
    if (aborted) {
      logInfo(seissol::MPI::mpi.rank()) << "Aborting simulation.";
      break;
    }

    // update the DOFs
    computeStopwatch.start();
    seissolInstance.timeManager().advanceInTime(upcomingTime);
    computeStopwatch.pause();

    ioStopwatch.start();

    // update current time
    currentTime = upcomingTime;

    // Synchronize data (TODO(David): synchronize lazily)
    seissolInstance.timeManager().synchronizeTo(seissol::initializer::AllocationPlace::Host);

    // Check all synchronization point hooks and set the new upcoming time
    upcomingTime = std::min(finalTime, Modules::callSyncHook(currentTime, timeTolerance));

    ioStopwatch.pause();

    double currentSplit = simulationStopwatch.split();
    Stopwatch::print("Time spent this phase (total):", currentSplit - lastSplit, seissol::MPI::mpi.comm());
    Stopwatch::print("Time spent this phase (compute):", computeStopwatch.split(), seissol::MPI::mpi.comm());
    Stopwatch::print("Time spent this phase (blocking IO):", ioStopwatch.split(), seissol::MPI::mpi.comm());
    seissolInstance.flopCounter().printPerformanceUpdate(currentSplit);
    lastSplit = currentSplit;
  }

  // synchronize data (TODO(David): synchronize lazily)
  seissolInstance.timeManager().synchronizeTo(seissol::initializer::AllocationPlace::Host);

  Modules::callSyncHook(currentTime, timeTolerance, true);

  const double wallTime = simulationStopwatch.pause();
  simulationStopwatch.printTime("Simulation time (total):", seissol::MPI::mpi.comm());
  computeStopwatch.printTime("Simulation time (compute):", seissol::MPI::mpi.comm());
  ioStopwatch.printTime("Simulation time (blocking IO):", seissol::MPI::mpi.comm());

  Modules::callHook<ModuleHook::SimulationEnd>();

  const auto& memoryManager = seissolInstance.getMemoryManager();
  const bool isLoopStatisticsNetcdfOutputOn = memoryManager.isLoopStatisticsNetcdfOutputOn();
  const auto& outputPrefix = memoryManager.getOutputPrefix();
  seissolInstance.timeManager().printComputationTime(outputPrefix, isLoopStatisticsNetcdfOutputOn);

  seissolInstance.analysisWriter().printAnalysis(currentTime);

  seissolInstance.flopCounter().printPerformanceSummary(wallTime);
}
