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

#include "Parallel/MPI.h"

#include <algorithm>
#include <cassert>
#include <cstring>

#include "FaultWriter.h"
#include "AsyncCellIDs.h"
#include "SeisSol.h"
#include "Modules/Modules.h"
#include "Solver/Interoperability.h"

extern seissol::Interoperability e_interoperability;

void seissol::writer::FaultWriter::setUp()
{
  setExecutor(m_executor);

  if (isAffinityNecessary()) {
    const auto freeCpus = SeisSol::main.getPinning().getFreeCPUsMask();
    logInfo(seissol::MPI::mpi.rank()) << "Fault writer thread affinity:"
      << parallel::Pinning::maskToString(freeCpus);
    if (parallel::Pinning::freeCPUsMaskEmpty(freeCpus)) {
      logError() << "There are no free CPUs left. Make sure to leave one for the I/O thread(s).";
    }
    setAffinityIfNecessary(freeCpus);
  }
}

void seissol::writer::FaultWriter::init(const unsigned int* cells, const double* vertices,
	unsigned int nCells, unsigned int nVertices,
	int* outputMask, const real** dataBuffer,
	const char* outputPrefix,
	double interval,
  xdmfwriter::BackendType backend)
{
	const int rank = seissol::MPI::mpi.rank();

	logInfo(rank) << "Initializing XDMF fault output.";

	// Initialize the asynchronous module
	async::Module<FaultWriterExecutor, FaultInitParam, FaultParam>::init();

	m_enabled = true;

	FaultInitParam param;
	param.timestep = m_timestep;
  param.backend = backend;

	// Create buffer for output prefix
	unsigned int bufferId = addSyncBuffer(outputPrefix, strlen(outputPrefix)+1, true);
	assert(bufferId == FaultWriterExecutor::OUTPUT_PREFIX); NDBG_UNUSED(bufferId);

	AsyncCellIDs<3> cellIds(nCells, nVertices, cells);

	// Create mesh buffers
	bufferId = addSyncBuffer(cellIds.cells(), nCells * 3 * sizeof(int));
	assert(bufferId == FaultWriterExecutor::CELLS);
	bufferId = addSyncBuffer(vertices, nVertices * 3 * sizeof(double));
	assert(bufferId == FaultWriterExecutor::VERTICES);

	// Create data buffers
	std::fill_n(param.outputMask, FaultInitParam::OUTPUT_MASK_SIZE, false);
	if (outputMask[0]) {
		param.outputMask[0] = true;
		param.outputMask[1] = true;
	}
	if (outputMask[1]) {
		param.outputMask[2] = true;
		param.outputMask[3] = true;
		param.outputMask[4] = true;
	}
	if (outputMask[2])
		param.outputMask[5] = true;
	if (outputMask[3]) {
		param.outputMask[6] = true;
		param.outputMask[7] = true;
	}
	if (outputMask[4]) {
		param.outputMask[8] = true;
		param.outputMask[9] = true;
		param.outputMask[10] = true;
	}
	if (outputMask[5]) {
		param.outputMask[11] = true;
		param.outputMask[12] = true;
	}
	if (outputMask[6])
		param.outputMask[13] = true;
	if (outputMask[7])
		param.outputMask[14] = true;
	if (outputMask[8])
		param.outputMask[15] = true;
	if (outputMask[9])
		param.outputMask[16] = true;
	if (outputMask[10])
		param.outputMask[17] = true;
        if (outputMask[11]) {
                param.outputMask[18] = true;
                param.outputMask[19] = true;
        }
	for (unsigned int i = 0; i < FaultInitParam::OUTPUT_MASK_SIZE; i++) {
		if (param.outputMask[i]) {
			addBuffer(dataBuffer[m_numVariables++], nCells * sizeof(real));
		}
	}

	//
	// Send all buffers for initialization
	//
	sendBuffer(FaultWriterExecutor::OUTPUT_PREFIX);

	sendBuffer(FaultWriterExecutor::CELLS);
	sendBuffer(FaultWriterExecutor::VERTICES);

	// Initialize the executor
	callInit(param);

	// Remove unused buffers
	removeBuffer(FaultWriterExecutor::OUTPUT_PREFIX);
	removeBuffer(FaultWriterExecutor::CELLS);
	removeBuffer(FaultWriterExecutor::VERTICES);

	// Register for the synchronization point hook
	Modules::registerHook(*this, SIMULATION_START);
	Modules::registerHook(*this, SYNCHRONIZATION_POINT);
	setSyncInterval(interval);
}

void seissol::writer::FaultWriter::simulationStart()
{
	syncPoint(0.0);
}

void seissol::writer::FaultWriter::syncPoint(double currentTime)
{
	SCOREP_USER_REGION("faultoutput_elementwise", SCOREP_USER_REGION_TYPE_FUNCTION)

	e_interoperability.calcElementwiseFaultoutput(currentTime);
	write(currentTime);
}
