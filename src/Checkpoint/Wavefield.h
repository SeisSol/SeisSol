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
 * Common interface for wave field checkpoints
 */

#ifndef CHECKPOINT_WAVEFIELD_H
#define CHECKPOINT_WAVEFIELD_H

#include "Parallel/MPI.h"

#include <cassert>

#include "utils/env.h"
#include "utils/logger.h"

#include "CheckPoint.h"
#include "WavefieldHeader.h"

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
	/** The header data (only used for loading) */
	WavefieldHeader* m_header;

	/** Pointer to the degrees of freedom */
	const real* m_dofs;

	/** Number of dofs */
	unsigned long m_numDofs;

	/** Number of (local) iterations we need to save all data (due to the 2GB limit) */
	unsigned int m_iterations;

	/** Number of total iteration we need to save the data (due to the 2GB limit) */
	unsigned int m_totalIterations;

	/** Number of cells that can be saved in one iteration (due to the 2GB limit) */
	const unsigned int m_dofsPerIteration;

public:
	Wavefield(unsigned long identifier)
		: CheckPoint(identifier),
		  m_header(0L),
		  m_dofs(0L), m_numDofs(0),
		  m_iterations(0), m_totalIterations(0),
		  m_dofsPerIteration((1ul<<30) / sizeof(real))
	{}

	virtual ~Wavefield() {}

	/**
	 * Set the header for loading checkpoints.
	 *
	 * This function can be overwritten to add custom header fields for the back-end
	 *
	 * @warning This function is only call on the compute nodes (not on I/O nodes)
	 * @todo One should probably split the checkpoint reader and writer
	 */
	virtual void setHeader(WavefieldHeader &header)
	{
		m_header = &header;
	}

	/**
	 * Initialize checkpointing
	 *
	 * @param groupSize The group size of the MPI group for asynchronous MPI checkpointing
	 *  (should be 1 for other modes)
	 *
	 * @return True of a valid checkpoint is available
	 */
	virtual bool init(size_t headerSize, unsigned long numDofs, unsigned int groupSize = 1)
	{
#ifndef USE_ASYNC_MPI
		assert(groupSize == 1);
#endif // USE_ASYNC_MPI

#ifdef USE_MPI
		// Setup rank, partitions, ...
		setComm(seissol::MPI::mpi.comm());
#endif // USE_MPI

		logInfo(rank()) << "Initializing check pointing";

		// Save size
		m_numDofs = numDofs;

		// Compute the group size/offset
		setGroupSumOffset(numDofs, groupSize);

		// Actual number of elements in the file (for this rank)
		unsigned long numDofsFile = numDofs;

#ifdef USE_MPI
		// More sure the local part is a multiple of the block size
		// This is important for all but the last rank
		int blockSize = utils::Env::get<int>("SEISSOL_CHECKPOINT_BLOCK_SIZE", 1);
		if (blockSize > 1 && rank() != partitions()-1) {
			if (blockSize % sizeof(real) != 0)
				logError() << "The block size for checkpointing must be a multiple of the size of the basic data type.";

			if (rank() % groupSize == groupSize - 1) {
				// Only the last rank of a group extends the number of elements
				unsigned int dofsPerBlock = blockSize / sizeof(real);
				unsigned long numGroupBlocks = (numGroupElems() + dofsPerBlock - 1) / dofsPerBlock;
				unsigned long numGroupDofsFile = numGroupBlocks * dofsPerBlock;
				numDofsFile = numDofs + (numGroupDofsFile - numGroupElems());
			}
		}
#endif // USE_MPI

		// Compute total number of cells and local offset
		setSumOffset(numDofsFile);

		// Work around 2 GB limit in MPI-IO
		m_iterations = (numDofs + m_dofsPerIteration - 1) / m_dofsPerIteration;
		m_totalIterations = m_iterations;
#ifdef USE_MPI
		MPI_Allreduce(MPI_IN_PLACE, &m_totalIterations, 1, MPI_UNSIGNED, MPI_MAX, comm());
#endif // USE_MPI

		return false;
	}

	/**
	 * Load a checkpoint file. Should only be done if init() returned true.
	 *
	 * @param dofs Pointer to the buffer where the degrees of freedom should be stored
	 */
	virtual void load(real* dofs) = 0;

	/**
	 * Initialize all const header data
	 */
	virtual void initHeader(WavefieldHeader &header)
	{
		header.identifier() = identifier();
	}

	/**
	 * Create checkpoint files. Should be called after loading the old checkpoint
	 */
	virtual void initLate(const real* dofs)
	{
		// Save dof pointer
		m_dofs = dofs;

		// Create checkpoint files
		createFiles();
	}

	/**
	 * Prepare writing a checkpoint
	 *
	 * Overwrite this function to start an asynchronous checkpoint write
	 *
	 * @param headerSize Size of the header in bytes
	 */
	virtual void writePrepare(const void* header, size_t headerSize)
	{
	}

	/**
	 * Write a checkpoint for the current time
	 *
	 * @param headerSize Size of the header in bytes
	 */
	virtual void write(const void* header, size_t headerSize) = 0;

protected:
	WavefieldHeader& header()
	{
		return *m_header;
	}

	bool hasHeader() const
	{
		return m_header != 0L;
	}

	void setDofs(const real* dofs)
	{
		m_dofs = dofs;
	}

	const real* dofs() const
	{
		return m_dofs;
	}

	unsigned long numDofs() const
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
