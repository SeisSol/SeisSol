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

#include "utils/env.h"
#include "utils/logger.h"

#include "Manager.h"
#include "SeisSol.h"

bool seissol::checkpoint::Manager::init(real* dofs, unsigned int numDofs,
		double* mu, double* slipRate1, double* slipRate2, double* slip, double* slip1, double* slip2,
		double* state, double* strength, unsigned int numSides, unsigned int numBndGP,
		int &faultTimeStep)
{
		if (m_backend == DISABLED) {
			// Always allocate the header struct because other still use it
			m_header.alloc();
			m_header.clear();
			return false;
		}

		// Initialize the asynchronous module
		async::Module<ManagerExecutor, CheckpointInitParam, CheckpointParam>::init();

		Wavefield* waveField;
		Fault* fault;
		createBackend(m_backend, waveField, fault);

		// Set the header
		waveField->setHeader(m_header);

		// Buffer for file name
		unsigned int id = addSyncBuffer(m_filename.c_str(), m_filename.size()+1, true);
		assert(id == FILENAME); NDBG_UNUSED(id);

		// Buffer for the header
		m_header.alloc(utils::Env::get<size_t>("SEISSOL_CHECKPOINT_ALIGNMENT", 0));
		id = addBuffer(m_header.data(), m_header.size(), true);
		assert(id == HEADER);

		// Buffers for data
		m_numDofs = numDofs;
		m_numDRDofs = numSides * numBndGP;

		id = addBuffer(dofs, numDofs * sizeof(real));
		assert(id == DOFS);
		id = addBuffer(mu, m_numDRDofs * sizeof(double));
		assert(id == DR_DOFS0);
		addBuffer(slipRate1, m_numDRDofs * sizeof(double));
		addBuffer(slipRate2, m_numDRDofs * sizeof(double));
		addBuffer(slip, m_numDRDofs * sizeof(double));
		addBuffer(slip1, m_numDRDofs * sizeof(double));
		addBuffer(slip2, m_numDRDofs * sizeof(double));
		addBuffer(state, m_numDRDofs * sizeof(double));
		addBuffer(strength, m_numDRDofs * sizeof(double));

		//
		// Initialization for loading checkpoints
		//
		waveField->setFilename(m_filename.c_str());
		fault->setFilename(m_filename.c_str());

		int exists = waveField->init(m_header.size(), numDofs, seissol::SeisSol::main.asyncIO().groupSize());
		exists &= fault->init(numSides, numBndGP,
			seissol::SeisSol::main.asyncIO().groupSize());

		// Make sure all ranks think the same about the existing checkpoint
#ifdef USE_MPI
		MPI_Allreduce(MPI_IN_PLACE, &exists, 1, MPI_INT, MPI_LAND, seissol::MPI::mpi.comm());
#endif // USE_MPI

		// Load checkpoint?
		if (exists) {
			waveField->load(dofs);
			fault->load(faultTimeStep, mu, slipRate1, slipRate2,
				slip, slip1, slip2, state, strength);
		} else {
			// Initialize header information (if not set from checkpoint)
			m_header.clear();
			waveField->initHeader(m_header);
		}

		waveField->close();
		fault->close();

		delete waveField;
		delete fault;

		sendBuffer(FILENAME,  m_filename.size()+1);

		// Initialize the executor
		CheckpointInitParam param;
		param.backend = m_backend;
		param.numBndGP = numBndGP;
		param.loaded = exists;
		callInit(param);

		removeBuffer(FILENAME);

		return exists;
}