// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

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

  // use an empty log message as visual separator
  logInfo() << "";

  while( m_finalTime > m_currentTime + timeTolerance ) {
    if (upcomingTime < m_currentTime + timeTolerance) {
      logError() << "Simulator did not advance in time from" << m_currentTime << "to" << upcomingTime;
    }
    if (m_abort) {
        logInfo() << "Aborting simulation.";
        break; 
    }

    // update the DOFs
    logInfo() << "Start simulation epoch.";
    computeStopwatch.start();
    seissolInstance.timeManager().advanceInTime( upcomingTime );
    computeStopwatch.pause();
    logInfo() << "End simulation epoch. Sync point.";

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

    // use an empty log message as visual separator
    logInfo() << "";
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

