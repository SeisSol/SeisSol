/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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

#ifndef FAULTWRITER_H
#define FAULTWRITER_H

#include "Parallel/MPI.h"

#include "utils/logger.h"

#include "async/Module.h"

#include "FaultWriterExecutor.h"
#include "Modules/Module.h"
#include "Monitoring/instrumentation.fpp"

namespace seissol
{

namespace writer
{

class FaultWriter : private async::Module<FaultWriterExecutor, FaultInitParam, FaultParam>,
	public seissol::Module
{
private:
	/** Is enabled? */
	bool m_enabled;

	/** The asynchronous executor */
	FaultWriterExecutor m_executor;

	/** Total number of variables */
	unsigned int m_numVariables;

public:
	FaultWriter()
		: m_enabled(false)
	{
	}

	/**
	 * Called by ASYNC on all ranks
	 */
	void setUp()
	{
		setExecutor(m_executor);
	}

	void init(const int* cells, const double* vertices,
		unsigned int nCells, unsigned int nVertices,
		int* outputMask, const double** dataBuffer,
		const char* outputPrefix,
		double interval);

	void write(double time)
	{
		SCOREP_USER_REGION("FaultWriter_write", SCOREP_USER_REGION_TYPE_FUNCTION)

		if (!m_enabled)
			logError() << "Trying to write fault output, but fault output is not enabled";

		const int rank = seissol::MPI::mpi.rank();

		wait();

		logInfo(rank) << "Writing faultoutput at time" << utils::nospace << time << ".";

		FaultParam param;
		param.time = time;

		for (unsigned int i = 0; i < m_numVariables; i++)
			sendBuffer(FaultWriterExecutor::VARIABLES0 + i);

		call(param);

		logInfo(rank) << "Writing faultoutput at time" << utils::nospace << time << ". Done.";
	}

	void close()
	{
		if (!m_enabled)
			return;

		wait();
	}

	void tearDown()
	{
		if (!m_enabled)
			return;

		m_executor.finalize();
	}

	//
	// Hooks
	//
	void simulationStart();

	void syncPoint(double currentTime);
};

}

}

#endif // FAULTWRITER_H