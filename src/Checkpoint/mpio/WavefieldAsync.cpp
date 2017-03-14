/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2016, SeisSol Group
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

#include <cstring>

#include "WavefieldAsync.h"

bool seissol::checkpoint::mpio::WavefieldAsync::init(unsigned long numDofs, unsigned int groupSize)
{
	bool exists = Wavefield::init(numDofs, groupSize);

	m_dofsCopy = new real[numDofs];

	return exists;
}

void seissol::checkpoint::mpio::WavefieldAsync::writePrepare(const void* header, size_t headerSize)
{
	EPIK_TRACER("CheckPoint_writePrepare");
	SCOREP_USER_REGION("CheckPoint_writePrepare", SCOREP_USER_REGION_TYPE_FUNCTION);

	// Write the header
	writeHeader(header, headerSize);

	// Create copy of the dofs
	memcpy(m_dofsCopy, dofs(), numDofs()*sizeof(real));

	// Save data
	EPIK_USER_REG(r_write_wavefield, "checkpoint_write_begin_wavefield");
	SCOREP_USER_REGION_DEFINE(r_write_wavefield);
	EPIK_USER_START(r_write_wavefield);
	SCOREP_USER_REGION_BEGIN(r_write_wavefield, "checkpoint_write_begin_wavefield", SCOREP_USER_REGION_TYPE_COMMON);

	checkMPIErr(setDataView(file()));

	checkMPIErr(MPI_File_write_all_begin(file(), m_dofsCopy, numDofs(), MPI_DOUBLE));

	EPIK_USER_END(r_write_wavefield);
	SCOREP_USER_REGION_END(r_write_wavefield);

	m_started = true;

	logInfo(rank()) << "Checkpoint backend: Writing. Done.";
}

void seissol::checkpoint::mpio::WavefieldAsync::write(const void* header, size_t headerSize)
{
	EPIK_TRACER("CheckPoint_write");
	SCOREP_USER_REGION("CheckPoint_write", SCOREP_USER_REGION_TYPE_FUNCTION);

	logInfo(rank()) << "Checkpoint backend: Writing.";

	if (m_started)
		checkMPIErr(MPI_File_write_all_end(file(), m_dofsCopy, MPI_STATUS_IGNORE));

	// Finalize the checkpoint
	finalizeCheckpoint();
}

void seissol::checkpoint::mpio::WavefieldAsync::close()
{
	// Finalize last checkpoint
	write(0, 0); // Time does not matter

	delete [] m_dofsCopy;

	Wavefield::close();
}
