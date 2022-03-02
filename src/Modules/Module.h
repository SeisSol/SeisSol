/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, SeisSol Group
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
 */

#ifndef MODULE_H
#define MODULE_H

#include <cassert>
#include <climits>
#include <cmath>
#include <limits>

#include "utils/logger.h"
#include "Parallel/MPI.h"

namespace seissol
{

/**
 * Base class for all modules
 */
class Module
{
private:
	/** The synchronization interval for this module */
	double m_syncInterval;

	/** The next synchronization point for this module */
	double m_nextSyncPoint;

  /** The last time when syncPoint was called */
	double m_lastSyncPoint;

public:
	/**
	 * Possible priorites for modules
	 *
	 * @todo Use std::numeric_limits<int> as soon as we switch to
	 *  C++0x
	 */
	enum Priority
	{
		MAX = INT_MIN,
		HIGHEST = -10000,
		HIGHER = -1000,
		HIGH = -100,
		DEFAULT = 0,
		LOW = 100,
		LOWER = 1000,
		LOWEST = 10000,
		MIN = INT_MAX
	};

public:
	Module()
		: m_syncInterval(0), m_nextSyncPoint(0), m_lastSyncPoint(-std::numeric_limits<double>::infinity())
	{ }

	/**
	 * Called by {@link Modules} at every synchronization point
	 *
	 * We have to ensure that this is "our" synchronization point before
	 * calling {@link syncPoint}.
	 *
	 * @return The next synchronization point for this module
	 */
	double potentialSyncPoint(double currentTime, double timeTolerance, bool forceSyncPoint)
	{
    if (std::abs(currentTime - m_lastSyncPoint) < timeTolerance) {
      int const rank = seissol::MPI::mpi.rank();
      logInfo(rank) << "Ignoring duplicate synchronisation point at time" << currentTime << "; the last sync point was at " << m_lastSyncPoint;
    } else if (forceSyncPoint || std::abs(currentTime - m_nextSyncPoint) < timeTolerance) {
			syncPoint(currentTime);
      m_lastSyncPoint = currentTime;
			m_nextSyncPoint += m_syncInterval;
		}

		return m_nextSyncPoint;
	}

	/**
	 * Called by {@link Modules} before the simulation starts to set the synchronization point.
	 *
	 * This is only called for modules that register for the SYNCHRONIZATION_POINT hook.
	 */
	void setSimulationStartTime(double time)
	{
		assert(m_syncInterval > 0);
		m_lastSyncPoint = time;
		m_nextSyncPoint = time + m_syncInterval;
	}

	//
	// Potential hooks
	//

	/**
	 * Called before initializing MPI
	 */
	virtual void preMPI()
	{
	}

	/**
	 * Called after MPI initialization
	 */
	virtual void postMPIInit()
	{
	}

	/**
	 * Called after mesh initialization
	 */
	virtual void postMesh()
	{
	}

	/**
	 * Called before the model is initialized
	 */
	virtual void preModel()
	{
	}

	/**
	 * Called after the model is initialized
	 */
	virtual void postModel()
	{
	}

	/**
	 * Called before the actual simulation.
	 *
	 * Only called when the simulation is not started from a checkpoint
	 */
	virtual void simulationStart()
	{
	}

	/**
	 * Called at synchronization points
	 */
	virtual void syncPoint(double currentTime)
	{
	}

protected:
  double syncInterval() const {
    return m_syncInterval;
  }

	/**
	 * Set the synchronization interval for this module
	 *
	 * This is only required for modules that register for {@link SYNCHRONIZATION_POINT}.
	 */
	void setSyncInterval(double interval)
	{
		if (m_syncInterval != 0)
			logError() << "Synchronization interval is already set";
		m_syncInterval = interval;
	}
};

}

#endif // MODULE_H
