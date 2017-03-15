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

#ifndef CHECKPOINT_MANAGER_EXECUTOR_H
#define CHECKPOINT_MANAGER_EXECUTOR_H

#include "async/ExecInfo.h"

#include "Backend.h"
#include "Monitoring/Stopwatch.h"

namespace seissol
{

namespace checkpoint
{

/** Buffer ids for asynchronous IO */
enum BufferTags {
	FILENAME = 0,
	HEADER = 1,
	DOFS = 2,
	DR_DOFS0 = 3
};

/**
 * Initialization parameters for checkpoints
 */
struct CheckpointInitParam
{
	Backend backend;
	unsigned int numBndGP;
	bool loaded;
};

/**
 * Parameters for checkpoints
 */
struct CheckpointParam
{
	double time;
	int faultTimeStep;
};

class ManagerExecutor
{
private:
	/** The wave field checkpoint */
	Wavefield *m_waveField;

	/** The dynamic rupture checkpoint */
	Fault *m_fault;

	/** Stopwatch for checkpoint backend */
	Stopwatch m_stopwatch;

public:
	ManagerExecutor()
		: m_waveField(0L),
		  m_fault(0L)
	{ }

	virtual ~ManagerExecutor()
	{ }

	/**
	 * Initializer on the executor
	 */
	void execInit(const async::ExecInfo &info, const CheckpointInitParam &param)
	{
		const char* filename = static_cast<const char*>(info.buffer(FILENAME));

		createBackend(param.backend, m_waveField, m_fault);

		m_waveField->setFilename(filename);
		m_fault->setFilename(filename);

		m_waveField->init(info.bufferSize(HEADER), info.bufferSize(DOFS) / sizeof(real));
		m_fault->init(info.bufferSize(DR_DOFS0) / param.numBndGP / sizeof(double), param.numBndGP);

		if (param.loaded) {
			m_waveField->setLoaded();
			m_fault->setLoaded();
		}

		const real* dofs = static_cast<const real*>(info.buffer(DOFS));
		const double* drDofs[8];
		for (unsigned int i = 0; i < 8; i++)
			drDofs[i] = static_cast<const double*>(info.buffer(DR_DOFS0 + i));

		m_waveField->initLate(dofs);
		m_fault->initLate(drDofs[0], drDofs[1], drDofs[2], drDofs[3], drDofs[4], drDofs[5],
			drDofs[6], drDofs[7]);
	}

	/**
	 * Write the checkpoint
	 */
	void exec(const async::ExecInfo &info, const CheckpointParam &param)
	{
		m_stopwatch.start();

		m_waveField->write(info.buffer(HEADER), info.bufferSize(HEADER));
		m_fault->write(param.faultTimeStep);

		// Update both links at the "same" time
		m_waveField->updateLink();
		m_fault->updateLink();

		// Prepare next checkpoint (only for async checkpoints)
		m_waveField->writePrepare(info.buffer(HEADER), info.bufferSize(HEADER));
		m_fault->writePrepare(param.faultTimeStep);

		m_stopwatch.pause();
	}

	void finalize()
	{
		if (m_waveField) {
			m_stopwatch.printTime("Time checkpoint backend:");

			m_waveField->close();
			m_fault->close();

			delete m_waveField;
			m_waveField = 0L;
			delete m_fault;
			m_fault = 0L;
		}
	}
};

}

}

#endif // CHECKPOINT_MANAGER_EXECUTOR_H