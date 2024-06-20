/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2017, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Entry point of the simulation.
 **/

#include <limits>

#include "Clustering/TimeManager.h"
#include "Modules/Modules.h"
#include "Monitoring/FlopCounter.hpp"
#include "Monitoring/Stopwatch.h"
#include "Monitoring/Unit.hpp"
#include "ResultWriter/AnalysisWriter.h"
#include "ResultWriter/EnergyOutput.h"
#include "SeisSol.h"
#include "Simulator.h"

seissol::Simulator::Simulator()
    : currentTime(0), finalTime(0), usePlasticity(false), checkPointTime(0),
      checkPointInterval(std::numeric_limits<double>::max()), aborted(false) {}

void seissol::Simulator::setCheckPointInterval(double checkPointInterval) {
  assert(checkPointInterval > 0);
  this->checkPointInterval = checkPointInterval;
}

bool seissol::Simulator::checkPointingEnabled() {
  return checkPointInterval < std::numeric_limits<double>::max();
}

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

  double timeTolerance = seissolInstance.timeManager().getTimeTolerance();

  // Write initial wave field snapshot
  if (currentTime == 0.0) {
    Modules::callHook<ModuleHook::SimulationStart>();
  }

  // intialize wave field and checkpoint time
  checkPointTime = currentTime;
  Modules::setSimulationStartTime(currentTime);

  // derive next synchronization time
  double upcomingTime = finalTime;
  // NOTE: This will not call the module specific implementation of the synchronization hook
  // since the current time is the simulation start time. We only use this function here to
  // get correct upcoming time. To be on the safe side, we use zero time tolerance.
  upcomingTime = std::min(upcomingTime, Modules::callSyncHook(currentTime, 0.0));
  upcomingTime = std::min(upcomingTime, std::abs(checkPointTime + checkPointInterval));

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

    // Set new upcoming time (might by overwritten by any of the modules)
    upcomingTime = finalTime;

    // synchronize data (TODO(David): synchronize lazily)
    seissolInstance.timeManager().synchronizeTo(seissol::initializer::AllocationPlace::Host);

    // Check all synchronization point hooks
    upcomingTime = std::min(upcomingTime, Modules::callSyncHook(currentTime, timeTolerance));

    // write checkpoint if required
    if (std::abs(currentTime - (checkPointTime + checkPointInterval)) < timeTolerance) {
      const unsigned int faultTimeStep = seissolInstance.faultWriter().timestep();
      seissolInstance.checkPointManager().write(currentTime, faultTimeStep);
      checkPointTime += checkPointInterval;
    }
    upcomingTime = std::min(upcomingTime, checkPointTime + checkPointInterval);

    ioStopwatch.pause();

    double currentSplit = simulationStopwatch.split();
    Stopwatch::print(
        "Time spent this phase (total):", currentSplit - lastSplit, seissol::MPI::mpi.comm());
    Stopwatch::print(
        "Time spent this phase (compute):", computeStopwatch.split(), seissol::MPI::mpi.comm());
    Stopwatch::print("Time spent this phase (IO):", ioStopwatch.split(), seissol::MPI::mpi.comm());
    seissolInstance.flopCounter().printPerformanceUpdate(currentSplit);
    lastSplit = currentSplit;
  }

  // synchronize data (TODO(David): synchronize lazily)
  seissolInstance.timeManager().synchronizeTo(seissol::initializer::AllocationPlace::Host);

  Modules::callSyncHook(currentTime, timeTolerance, true);

  double wallTime = simulationStopwatch.pause();
  simulationStopwatch.printTime("Simulation time (total):", seissol::MPI::mpi.comm());
  computeStopwatch.printTime("Simulation time (compute):", seissol::MPI::mpi.comm());
  ioStopwatch.printTime("Simulation time (IO):", seissol::MPI::mpi.comm());

  Modules::callHook<ModuleHook::SimulationEnd>();

  const auto& memoryManager = seissolInstance.getMemoryManager();
  const bool isLoopStatisticsNetcdfOutputOn = memoryManager.isLoopStatisticsNetcdfOutputOn();
  const auto& outputPrefix = memoryManager.getOutputPrefix();
  seissolInstance.timeManager().printComputationTime(outputPrefix, isLoopStatisticsNetcdfOutputOn);

  seissolInstance.analysisWriter().printAnalysis(currentTime);

  seissolInstance.flopCounter().printPerformanceSummary(wallTime);
}
