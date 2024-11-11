/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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

#include "Simulator.h"
#include "Modules/Modules.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Stopwatch.h"
#include "Monitoring/Unit.h"
#include "ResultWriter/AnalysisWriter.h"
#include "ResultWriter/EnergyOutput.h"
#include "SeisSol.h"
#include "time_stepping/TimeManager.h"

seissol::Simulator::Simulator():
  m_currentTime(        0 ),
  m_finalTime(          0 ),
  m_usePlasticity(  false ),
  m_abort( false ) {}

void seissol::Simulator::setFinalTime( double i_finalTime ) {
  assert( i_finalTime > 0 );
  m_finalTime = i_finalTime;
}

void seissol::Simulator::setUsePlasticity( bool plasticity ) {
  m_usePlasticity = plasticity;
}

void seissol::Simulator::setCurrentTime( double i_currentTime ) {
	assert( i_currentTime >= 0 );
	m_currentTime = i_currentTime;
}

void seissol::Simulator::abort() {
	m_abort = true;
}


void seissol::Simulator::simulate(seissol::SeisSol& seissolInstance) {
  SCOREP_USER_REGION( "simulate", SCOREP_USER_REGION_TYPE_FUNCTION )

  auto* faultOutputManager = seissolInstance.timeManager().getFaultOutputManager();
  faultOutputManager->writePickpointOutput(0.0, 0.0);

  Stopwatch simulationStopwatch;
  simulationStopwatch.start();

  Stopwatch computeStopwatch;
  Stopwatch ioStopwatch;

  ioStopwatch.start();

  // Set start time (required for checkpointing)
  seissolInstance.timeManager().setInitialTimes(m_currentTime);

  double timeTolerance = seissolInstance.timeManager().getTimeTolerance();

  // Write initial wave field snapshot
  if (m_currentTime == 0.0) {
    Modules::callHook<ModuleHook::SimulationStart>();
  }

  // intialize wave field and checkpoint time
  Modules::setSimulationStartTime(m_currentTime);

  // derive next synchronization time
  double upcomingTime = m_finalTime;
  // NOTE: This will not call the module specific implementation of the synchronization hook
  // since the current time is the simulation start time. We only use this function here to
  // get correct upcoming time. To be on the safe side, we use zero time tolerance.
  upcomingTime = std::min( upcomingTime, Modules::callSyncHook(m_currentTime, 0.0) );

  double lastSplit = 0;

  ioStopwatch.pause();

  Stopwatch::print("Time spent for initial IO:", ioStopwatch.split(), seissol::MPI::mpi.comm());

  while( m_finalTime > m_currentTime + timeTolerance ) {
    if (upcomingTime < m_currentTime + timeTolerance) {
      logError() << "Simulator did not advance in time from" << m_currentTime << "to" << upcomingTime;
    }
    if (m_abort) {
        logInfo() << "Aborting simulation.";
        break; 
    }

    // update the DOFs
    computeStopwatch.start();
    seissolInstance.timeManager().advanceInTime( upcomingTime );
    computeStopwatch.pause();

    ioStopwatch.start();

    // update current time
    m_currentTime = upcomingTime;

    // Synchronize data (TODO(David): synchronize lazily)
    seissolInstance.timeManager().synchronizeTo(seissol::initializer::AllocationPlace::Host);

    // Check all synchronization point hooks and set the new upcoming time
    upcomingTime = std::min(m_finalTime, Modules::callSyncHook(m_currentTime, timeTolerance));

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

  Modules::callSyncHook(m_currentTime, timeTolerance, true);

  double wallTime = simulationStopwatch.pause();
  simulationStopwatch.printTime("Simulation time (total):", seissol::MPI::mpi.comm());
  computeStopwatch.printTime("Simulation time (compute):", seissol::MPI::mpi.comm());
  ioStopwatch.printTime("Simulation time (blocking IO):", seissol::MPI::mpi.comm());

  Modules::callHook<ModuleHook::SimulationEnd>();

  const auto& memoryManager = seissolInstance.getMemoryManager();
  const bool isLoopStatisticsNetcdfOutputOn = memoryManager.isLoopStatisticsNetcdfOutputOn();
  const auto& outputPrefix = memoryManager.getOutputPrefix();
  seissolInstance.timeManager().printComputationTime(outputPrefix,
                                                            isLoopStatisticsNetcdfOutputOn);

  seissolInstance.analysisWriter().printAnalysis(m_currentTime);

  seissolInstance.flopCounter().printPerformanceSummary(wallTime);
}
