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
 */

#include <mpi.h>

#include "utils/env.h"

#include "Wavefield.h"
#include "Monitoring/instrumentation.fpp"

void seissol::checkpoint::mpio::Wavefield::setHeader(seissol::checkpoint::WavefieldHeader &header)
{
	seissol::checkpoint::Wavefield::setHeader(header);
	header.add(m_partitionComp);
}

bool seissol::checkpoint::mpio::Wavefield::init(size_t headerSize, unsigned long numDofs, unsigned int groupSize)
{
	seissol::checkpoint::Wavefield::init(headerSize, numDofs, groupSize);

	// Create the header data type
	// We cannot use header since this will be called on I/O nodes as well
	MPI_Datatype headerType;
	MPI_Type_contiguous(headerSize, MPI_BYTE, &headerType);
	setHeaderType(headerType);

	// Define the file view
	defineFileView(headerSize, sizeof(real), numDofs);

	// Large buffers are supported?
	m_useLargeBuffer = utils::Env::get<int>("SEISSOL_CHECKPOINT_MPIO_LARGE_BUFFER", 1) != 0;

	return exists();
}

void seissol::checkpoint::mpio::Wavefield::load(real* dofs)
{
	logInfo(rank()) << "Loading wave field checkpoint";

	seissol::checkpoint::CheckPoint::setLoaded();

	MPI_File file = open();
	if (file == MPI_FILE_NULL)
		logError() << "Could not open checkpoint file";

	// Read and broadcast header
	checkMPIErr(setHeaderView(file));

	if (rank() == 0)
		checkMPIErr(MPI_File_read(file, header().data(), 1, headerType(), MPI_STATUS_IGNORE));

	MPI_Bcast(header().data(), 1, headerType(), 0, comm());

	// Read dofs
	checkMPIErr(setDataView(file));
	checkMPIErr(MPI_File_read_all(file, dofs, numDofs(), MPI_DOUBLE, MPI_STATUS_IGNORE));

	// Close the file
	checkMPIErr(MPI_File_close(&file));
}

void seissol::checkpoint::mpio::Wavefield::initHeader(WavefieldHeader &header)
{
	seissol::checkpoint::Wavefield::initHeader(header);

	header.value(m_partitionComp) = partitions();
}

void seissol::checkpoint::mpio::Wavefield::write(const void* header, size_t headerSize)
{
	SCOREP_USER_REGION("CheckPoint_write", SCOREP_USER_REGION_TYPE_FUNCTION);

	logInfo(rank()) << "Checkpoint backend: Writing.";

	// Write the header
	writeHeader(header, headerSize);

	// Save data
	SCOREP_USER_REGION_DEFINE(r_write_wavefield);
	SCOREP_USER_REGION_BEGIN(r_write_wavefield, "checkpoint_write_wavefield", SCOREP_USER_REGION_TYPE_COMMON);
	checkMPIErr(setDataView(file()));

	unsigned int totalIter = totalIterations();
	unsigned int iter = iterations();
	unsigned int count = dofsPerIteration();
	if (m_useLargeBuffer) {
		totalIter = (totalIter + sizeof(real) - 1) / sizeof(real);
		iter = (iter + sizeof(real) - 1) / sizeof(real);
		count *= sizeof(real);
	}
	unsigned long offset = 0;
	for (unsigned int i = 0; i < totalIter; i++) {
		if (i == iter-1)
			// Last iteration
			count = numDofs() - (iter-1) * count;

		checkMPIErr(MPI_File_write_all(file(), const_cast<real*>(&dofs()[offset]), count, MPI_DOUBLE, MPI_STATUS_IGNORE));

		if (i < iter-1)
			offset += count;
		// otherwise we just continue writing the last chunk over and over
		else if (i != totalIter-1)
			checkMPIErr(MPI_File_seek(file(), -count * sizeof(real), MPI_SEEK_CUR));
	}

	SCOREP_USER_REGION_END(r_write_wavefield);

	// Finalize the checkpoint
	finalizeCheckpoint();

	logInfo(rank()) << "Checkpoint backend: Writing. Done.";
}

bool seissol::checkpoint::mpio::Wavefield::validate(MPI_File file)
{
	if (setHeaderView(file) != 0) {
		logWarning() << "Could not set checkpoint header view";
		return false;
	}

	int result = true;

	if (rank() == 0 && hasHeader()) { // Only validate on compute nodes
		// Check the header
		MPI_File_read(file, header().data(), 1, headerType(), MPI_STATUS_IGNORE);

		if (header().identifier() != identifier()) {
			logWarning() << "Checkpoint identifier does match";
			result = false;
		} else if (header().value(m_partitionComp) != partitions()) {
			logWarning() << "Number of partitions in checkpoint does not match";
			result = false;
		}
	}

	// Make sure everybody knows the result of the validation
	MPI_Bcast(&result, 1, MPI_INT, 0, comm());

	return result;
}

void seissol::checkpoint::mpio::Wavefield::writeHeader(const void* header, size_t headerSize)
{
	SCOREP_USER_REGION("checkpoint_write_header", SCOREP_USER_REGION_TYPE_FUNCTION);

	checkMPIErr(setHeaderView(file()));

	if (rank() == 0)
		checkMPIErr(MPI_File_write(file(), const_cast<void*>(header), 1, headerType(), MPI_STATUS_IGNORE));

}