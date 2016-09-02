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
 */

#include <mpi.h>

#include <cstring>

#include "FaultAsync.h"

bool seissol::checkpoint::mpio::FaultAsync::init(unsigned int numSides, unsigned int numBndGP,
		unsigned int groupSize)
{
	bool exists = Fault::init(numSides, numBndGP, groupSize);

	if (numSides != 0)
		m_dataCopy = new double[NUM_VARIABLES * numSides * numBndGP];

	return exists;
}

void seissol::checkpoint::mpio::FaultAsync::writePrepare(int timestepFault)
{
	EPIK_TRACER("CheckPointFault_writePrepare");
	SCOREP_USER_REGION("CheckPointFault_writePrepare", SCOREP_USER_REGION_TYPE_FUNCTION);

	if (numSides() == 0)
		return;

	// Write the header
	writeHeader(timestepFault);

	// Create copy of the data
	for (unsigned int i = 0; i < NUM_VARIABLES; i++)
		memcpy(&m_dataCopy[i*numSides()*numBndGP()],
				data(i), numSides()*numBndGP()*sizeof(double));

	// Save data
	EPIK_USER_REG(r_write_wavefield, "checkpoint_write_begin_fault");
	SCOREP_USER_REGION_DEFINE(r_write_fault);
	EPIK_USER_START(r_write_wavefield);
	SCOREP_USER_REGION_BEGIN(r_write_fault, "checkpoint_write_begin_fault", SCOREP_USER_REGION_TYPE_COMMON);

	checkMPIErr(setDataView(file()));
	checkMPIErr(MPI_File_write_all_begin(file(), m_dataCopy, numSides() * numBndGP() * NUM_VARIABLES, MPI_DOUBLE));

	EPIK_USER_END(r_write_fault);
	SCOREP_USER_REGION_END(r_write_fault);

	m_started = true;

	logInfo(rank()) << "Checkpoint backend: Writing fault. Done.";
}

void seissol::checkpoint::mpio::FaultAsync::write(int timestepFault)
{
	EPIK_TRACER("CheckPointFault_write");
	SCOREP_USER_REGION("CheckPointFault_write", SCOREP_USER_REGION_TYPE_FUNCTION);

	if (numSides() == 0)
		return;

	logInfo(rank()) << "Checkpoint backend: Writing fault.";

	if (m_started)
		checkMPIErr(MPI_File_write_all_end(file(), m_dataCopy, MPI_STATUS_IGNORE));

	// Finalize the checkpoint
	finalizeCheckpoint();
}

void seissol::checkpoint::mpio::FaultAsync::close()
{
	if (numSides() != 0) {
		// Finalize last checkpoint
		write(0); // Time does not matter

		delete [] m_dataCopy;
	}

	Fault::close();
}
