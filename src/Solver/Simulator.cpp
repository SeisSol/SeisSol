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
#include "SeisSol.h"
#include "Interoperability.h"
#include "time_stepping/TimeManager.h"
#include "Modules/Modules.h"
#include "Monitoring/Stopwatch.h"
#include "Monitoring/FlopCounter.hpp"
#include "ResultWriter/AnalysisWriter.h"
#include "ResultWriter/EnergyOutput.h"

extern seissol::Interoperability e_interoperability;

seissol::Simulator::Simulator():
  m_currentTime(        0 ),
  m_finalTime(          0 ),
  m_checkPointTime(     0 ),
  m_checkPointInterval( std::numeric_limits< double >::max() ),
  m_loadCheckPoint( false ) {}

void seissol::Simulator::setCheckPointInterval( double i_checkPointInterval ) {
  assert( m_checkPointInterval > 0 );
  m_checkPointInterval = i_checkPointInterval;
}

bool seissol::Simulator::checkPointingEnabled() {
  return m_checkPointInterval < std::numeric_limits<double>::max();
}

void seissol::Simulator::setFinalTime( double i_finalTime ) {
  assert( i_finalTime > 0 );
  m_finalTime = i_finalTime;
}

void seissol::Simulator::setCurrentTime( double i_currentTime ) {
	assert( i_currentTime > 0 );
	m_currentTime = i_currentTime;
}

void seissol::Simulator::simulate() {
  SCOREP_USER_REGION( "simulate", SCOREP_USER_REGION_TYPE_FUNCTION )

  Stopwatch stopwatch;
  stopwatch.start();

  // Set start time (required for checkpointing)
  seissol::SeisSol::main.timeManager().setInitialTimes(m_currentTime);

  // tolerance in time which is neglected
  double l_timeTolerance = seissol::SeisSol::main.timeManager().getTimeTolerance();

  // Copy initial dynamic rupture in order to ensure correct initial fault output
  e_interoperability.copyDynamicRuptureState();

  // Write initial wave field snapshot
  if (m_currentTime == 0.0) {
    Modules::callHook<SIMULATION_START>();
  }

  // intialize wave field and checkpoint time
  m_checkPointTime = m_currentTime;
  Modules::setSimulationStartTime(m_currentTime);

  // start the communication thread (if applicable)
  seissol::SeisSol::main.timeManager().startCommunicationThread();

  // derive next synchronization time
  double upcomingTime = m_finalTime;
  // NOTE: This will not call the module specific implementation of the synchronization hook
  // since the current time is the simulation start time. We only use this function here to
  // get correct upcoming time. To be on the safe side, we use zero time tolerance.
  upcomingTime = std::min( upcomingTime, Modules::callSyncHook(m_currentTime, 0.0) );
  upcomingTime = std::min( upcomingTime, std::abs(m_checkPointTime + m_checkPointInterval) );

  while( m_finalTime > m_currentTime + l_timeTolerance ) {
    if (upcomingTime < m_currentTime + l_timeTolerance)
      logError() << "Simulator did not advance in time from" << m_currentTime << "to" << upcomingTime;

    // update the DOFs
    seissol::SeisSol::main.timeManager().advanceInTime( upcomingTime );

    // update current time
    m_currentTime = upcomingTime;

    // Set new upcoming time (might by overwritten by any of the modules)
    upcomingTime = m_finalTime;

    // Check all synchronization point hooks
    upcomingTime = std::min(upcomingTime, Modules::callSyncHook(m_currentTime, l_timeTolerance));

    // write checkpoint if required
    if( std::abs( m_currentTime - ( m_checkPointTime + m_checkPointInterval ) ) < l_timeTolerance ) {
      const unsigned int faultTimeStep = seissol::SeisSol::main.faultWriter().timestep();
      seissol::SeisSol::main.checkPointManager().write(m_currentTime, faultTimeStep);
      m_checkPointTime += m_checkPointInterval;
    }
    upcomingTime = std::min(upcomingTime, m_checkPointTime + m_checkPointInterval);

    printPerformance(stopwatch.split());
  }
  
  Modules::callSyncHook(m_currentTime, l_timeTolerance, true);

  // stop the communication thread (if applicable)
  seissol::SeisSol::main.timeManager().stopCommunicationThread();

  double wallTime = stopwatch.split();
  logInfo(seissol::MPI::mpi.rank()) << "Elapsed time (via clock_gettime):" << wallTime << "seconds.";

  seissol::SeisSol::main.timeManager().printComputationTime();

  seissol::SeisSol::main.analysisWriter().printAnalysis(m_currentTime);

  printFlops();

  MeshReader& meshReader = seissol::SeisSol::main.meshReader();
  auto ltsTree = seissol::SeisSol::main.getMemoryManager().getLtsTree();
  auto lts     = seissol::SeisSol::main.getMemoryManager().getLts();
  auto* ltsLut = e_interoperability.getLtsLut();
  seissol::writer::printPlasticMoment(meshReader, ltsTree, lts, ltsLut);


}
