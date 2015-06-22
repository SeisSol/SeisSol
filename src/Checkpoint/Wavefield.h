/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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
 * Common interface for wave field checkpoints
 */

#ifndef CHECKPOINT_WAVEFIELD_H
#define CHECKPOINT_WAVEFIELD_H

#include "utils/env.h"
#include "utils/logger.h"

#include "CheckPoint.h"

#include "Initializer/preProcessorMacros.fpp"
#include "Initializer/typedefs.hpp"

namespace seissol
{

namespace checkpoint
{

/**
 * Common interface for wave field checkpoints
 */
class Wavefield : virtual public CheckPoint
{
private:
	/** Pointer to the degrees of freedom */
	real* m_dofs;

	/** Number of dofs */
	unsigned int m_numDofs;

	/** Number of (local) iterations we need to save all data (due to the 2GB limit) */
	unsigned int m_iterations;

	/** Number of total iteration we need to save the data (due to the 2GB limit) */
	unsigned int m_totalIterations;

	/** Number of cells that can be saved in one iteration (due to the 2GB limit) */
	const unsigned int m_dofsPerIteration;

public:
	Wavefield()
		: m_dofs(0L), m_numDofs(0),
		  m_iterations(0), m_totalIterations(0),
		  m_dofsPerIteration((1ul<<30) / sizeof(real))
	{}

	virtual ~Wavefield() {}

	/**
	 * Initialize checkpointing
	 *
	 * @return True of a valid checkpoint is available
	 */
	virtual bool init(real* dofs, unsigned int numDofs)
	{
#ifdef USE_MPI
		// Setup rank, partitions, ...
		setComm(MPI_COMM_WORLD);
#endif // USE_MPI

		logInfo(rank()) << "Initializing check pointing";

		// Save the dof pointer and size
		m_dofs = dofs;
		m_numDofs = numDofs;

		// Actual number of elements in the file (for this rank)
		unsigned int numDofsFile = numDofs;

#ifdef USE_MPI
		// More sure the local part is a multiple of the block size
		// This is important for all but the last rank
		unsigned int blockSize = utils::Env::get<unsigned int>("SEISSOL_CHECKPOINT_BLOCK_SIZE", 1);
		if (blockSize > 1 && rank() != partitions()-1) {
			if (blockSize % sizeof(real) != 0)
				logError() << "The block size for checkpointing must be a multiple of the size of the basic data type.";

			unsigned int dofsPerBlock = blockSize / sizeof(real);
			unsigned int numBlocks = (numDofs + dofsPerBlock - 1) / dofsPerBlock;
			numDofsFile = numBlocks * dofsPerBlock;
		}
#endif // USE_MPI

		// Compute total number of cells and local offset
		setSumOffset(numDofsFile);

		// Work around 2 GB limit in MPI-IO
		m_iterations = (numDofs + m_dofsPerIteration - 1) / m_dofsPerIteration;
		m_totalIterations = m_iterations;
#ifdef USE_MPI
		MPI_Allreduce(MPI_IN_PLACE, &m_totalIterations, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
#endif // USE_MPI

		return false;
	}

	/**
	 * Load a checkpoint file. Should only be done if init() returned true.
	 *
	 * @param[out] time Time of the simulation in the checkpoint
	 * @param[out] timestepWavefield Time step of the wave field writer in the checkpoint
	 *  (if the wave field writer was active)
	 */
	virtual void load(double &time, int &timestepWavefield) = 0;

	/**
	 * Write a checkpoint for the current time
	 *
	 * @param time The current time
	 * @param timestepWaveField The current time step of the wave field
	 */
	virtual void write(double time, int timestepWaveField) = 0;

protected:
	real* dofs()
	{
		return m_dofs;
	}

	const real* dofs() const
	{
		return m_dofs;
	}

	unsigned int numDofs() const
	{
		return m_numDofs;
	}

	unsigned int iterations() const
	{
		return m_iterations;
	}

	unsigned int totalIterations() const
	{
		return m_totalIterations;
	}

	unsigned int dofsPerIteration() const
	{
		return m_dofsPerIteration;
	}
};

}

}

#endif // CHECKPOINT_WAVEFIELD_H
