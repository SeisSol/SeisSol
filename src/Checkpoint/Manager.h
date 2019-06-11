/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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
 * Manages the checkpoints
 */

#ifndef CHECKPOINT_MANAGER_H
#define CHECKPOINT_MANAGER_H

#include "Parallel/MPI.h"
#include "Parallel/Pin.h"

#include <cassert>
#include <cstring>
#include <string>

#include "utils/logger.h"

#include "async/Module.h"

#include "Backend.h"
#include "ManagerExecutor.h"
#include "Wavefield.h"
#include "Fault.h"
#include "WavefieldHeader.h"
#include "Monitoring/Stopwatch.h"

namespace seissol
{

namespace checkpoint
{

class Manager : private async::Module<ManagerExecutor, CheckpointInitParam, CheckpointParam>
{
private:
	ManagerExecutor m_executor;

	/** The backend that should be used */
	Backend m_backend;

	/** The filename for the checkpoints */
	std::string m_filename;

	/** Number of DOFs */
	unsigned int m_numDofs;

	/** Number of DR DOFs */
	unsigned int m_numDRDofs;

	/** Checkpoint header */
	WavefieldHeader m_header;

	/** Stopwatch for checkpointing frontend */
	Stopwatch m_stopwatch;

public:
	Manager()
		: m_backend(DISABLED),
		  m_numDofs(0), m_numDRDofs(0)
	{
	}

	virtual ~Manager()
	{ }

	void setBackend(Backend backend)
	{
		m_backend = backend;
	}

	/**
	 * Set the filename prefix for checkpointing
	 *
	 * @param filename The filename prefix
	 */
	void setFilename(const char* filename)
	{
		m_filename = filename;
	}


	/**
	 * This is called on all ranks
	 */
	void setUp()
	{
		setExecutor(m_executor);
		if (isAffinityNecessary()) {
		  const auto freeCpus = parallel::getFreeCPUsMask();
		  logInfo(seissol::MPI::mpi.rank()) << "Checkpoint thread affinity:" << parallel::maskToString(parallel::getFreeCPUsMask());
		  if (parallel::freeCPUsMaskEmpty(freeCpus)) {
		    logError() << "There are no free CPUs left. Make sure to leave one for the I/O thread(s).";
		  }
		  setAffinityIfNecessary(freeCpus);
		}
	}

	/**
	 * The header struct
	 */
	WavefieldHeader& header()
	{
		return m_header;
	}

	/**
	 * Initialize checkpointing and load the last checkpoint if present
	 *
	 * @return True is a checkpoint was loaded, false otherwise
	 */
	bool init(real* dofs, unsigned int numDofs,
			double* mu, double* slipRate1, double* slipRate2, double* slip, double* slip1, double* slip2,
			double* state, double* strength, unsigned int numSides, unsigned int numBndGP,
			int &faultTimeStep);

	/**
	 * Write a checkpoint for the current time
	 *
	 * @param time The current time
	 * @param faultTimeStep The time step of the fault writer
	 */
	void write(double time, int faultTimeStep)
	{
		SCOREP_USER_REGION("CheckpointManager_write", SCOREP_USER_REGION_TYPE_FUNCTION);

		if (m_backend == DISABLED)
			return;

		m_stopwatch.start();

		const int rank = seissol::MPI::mpi.rank();

		// Set current time
		m_header.time() = time;

		SCOREP_USER_REGION_DEFINE(r_wait);
		SCOREP_USER_REGION_BEGIN(r_wait, "checkpointmanager_wait", SCOREP_USER_REGION_TYPE_COMMON);
		logInfo(rank) << "Checkpoint: Waiting for last.";
		wait();
		SCOREP_USER_REGION_END(r_wait);

		logInfo(rank) << "Checkpoint: Writing at time" << utils::nospace << time << '.';

		// Send buffers
		sendBuffer(HEADER);
		sendBuffer(DOFS, m_numDofs * sizeof(real));
		for (unsigned int i = 0; i < 8; i++)
			sendBuffer(DR_DOFS0+i, m_numDRDofs * sizeof(double));


		SCOREP_USER_REGION_DEFINE(r_call);
		SCOREP_USER_REGION_BEGIN(r_call, "checkpointmanager_call", SCOREP_USER_REGION_TYPE_COMMON);
		CheckpointParam param;
		param.time = time;
		param.faultTimeStep = faultTimeStep;
		call(param);
		SCOREP_USER_REGION_END(r_call);

		m_stopwatch.pause();

		logInfo(rank) << "Checkpoint: Writing at time" << utils::nospace << time << ". Done.";
	}

	/**
	 * Close checkpointing
	 */
	void close()
	{
		if (m_backend == DISABLED)
			return;

		// Terminate the executor
		wait();

		m_stopwatch.printTime("Time checkpoint frontend:");

		// Cleanup the asynchronous module
		async::Module<ManagerExecutor, CheckpointInitParam, CheckpointParam>::finalize();
	}

	/**
	 * Called on all ranks
	 */
	void tearDown()
	{
		m_executor.finalize();
	}

private:
};

}

}

#endif // CHECKPOINT_MANAGER_H
