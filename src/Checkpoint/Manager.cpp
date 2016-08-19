/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
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

#include "Manager.h"
#include "SeisSol.h"


bool seissol::checkpoint::Manager::init(real* dofs, unsigned int numDofs,
		double* mu, double* slipRate1, double* slipRate2, double* slip, double* slip1, double* slip2,
		double* state, double* strength, unsigned int numSides, unsigned int numBndGP,
		double &time, int &waveFieldTimeStep, int &faultTimeStep)
{
	// Set the executor
	setUp();

	if (m_backend == DISABLED)
		return false;

	// Buffer for file name
#ifdef USE_ASYNC_MPI
	m_bufferIds[FILENAME] = addBuffer(m_filename.size()+1);
#endif // USE_ASYNC_MPI

	// Buffers for data
#if defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
	m_bufferIds[DOFS] = addBuffer(numDofs * sizeof(real));
	m_bufferIds[MU] = addBuffer(numSides * numBndGP * sizeof(double));
	m_bufferIds[SLIP_RATE1] = addBuffer(numSides * numBndGP * sizeof(double));
	m_bufferIds[SLIP_RATE2] = addBuffer(numSides * numBndGP * sizeof(double));
	m_bufferIds[SLIP] = addBuffer(numSides * numBndGP * sizeof(double));
	m_bufferIds[SLIP1] = addBuffer(numSides * numBndGP * sizeof(double));
	m_bufferIds[SLIP2] = addBuffer(numSides * numBndGP * sizeof(double));
	m_bufferIds[STATE] = addBuffer(numSides * numBndGP * sizeof(double));
	m_bufferIds[STRENGTH] = addBuffer(numSides * numBndGP * sizeof(double));
#endif // defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)

	// Save pointers for later
	m_dofs = dofs;
	m_drDofs[0] = mu;
	m_drDofs[1] = slipRate1;
	m_drDofs[2] = slipRate2;
	m_drDofs[3] = slip;
	m_drDofs[4] = slip1;
	m_drDofs[5] = slip2;
	m_drDofs[6] = state;
	m_drDofs[7] = strength;

	m_numDofs = numDofs;
	m_numDRDofs = numSides * numBndGP;

	initBackend();

	m_waveField->setFilename(m_filename.c_str());
	m_fault->setFilename(m_filename.c_str());

	int exists = m_waveField->init(numDofs, seissol::SeisSol::main.asyncIO().groupSize());
	exists &= m_fault->init(numSides, numBndGP,
			seissol::SeisSol::main.asyncIO().groupSize());

	// Make sure all rank think the same about the existing checkpoint
#ifdef USE_MPI
	MPI_Allreduce(MPI_IN_PLACE, &exists, 1, MPI_INT, MPI_LAND, seissol::MPI::mpi.comm());
#endif // USE_MPI

	// Load checkpoint?
	if (exists) {
		m_waveField->load(time, waveFieldTimeStep, dofs);
		m_fault->load(faultTimeStep, mu, slipRate1, slipRate2,
			slip, slip1, slip2, state, strength);
	}

#ifdef USE_ASYNC_MPI
	fillBuffer(m_bufferIds[FILENAME], m_filename.c_str(), m_filename.size()+1);
#endif // USE_ASYNC_MPI

	// Initialize the executor
	CheckpointInitParam param;
#ifdef USE_ASYNC_MPI
	assert(sizeof(param.bufferIds) == sizeof(m_bufferIds));
	memcpy(param.bufferIds, m_bufferIds, sizeof(m_bufferIds));
	param.backend = m_backend;
#endif // USE_ASYNC_MPI
	param.numBndGP = numBndGP;
	param.loaded = exists;
	callInit(param);

	return exists;
}