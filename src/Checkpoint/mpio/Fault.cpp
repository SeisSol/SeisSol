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

#include <cstddef>

#include "Fault.h"

bool seissol::checkpoint::mpio::Fault::init(unsigned int numSides, unsigned int numBndGP,
		unsigned int groupSize)
{
	seissol::checkpoint::Fault::init(numSides, numBndGP, groupSize);

	if (numSides == 0)
		return true;

	// Compute total number of cells and local offset
	setSumOffset(numSides);

	// Create the header data type
	MPI_Datatype headerType;
	int blockLength[] = {1, 1};
	MPI_Aint displ[] = {offsetof(Header, identifier), offsetof(Header, timestepFault)};
	MPI_Datatype types[] = {MPI_UNSIGNED_LONG, MPI_INT};
	MPI_Type_create_struct(2, blockLength, displ, types, &headerType);
	setHeaderType(headerType);

	// Define the file view
	defineFileView(sizeof(Header), numBndGP * sizeof(double), numSides, NUM_VARIABLES);

	return exists();
}

void seissol::checkpoint::mpio::Fault::load(int &timestepFault, double* mu, double* slipRate1, double* slipRate2,
	double* slip, double* slip1, double* slip2, double* state, double* strength)
{
	if (numSides() == 0)
		return;

	logInfo(rank()) << "Loading fault checkpoint";

	seissol::checkpoint::CheckPoint::setLoaded();

	MPI_File file = open();
	if (file == MPI_FILE_NULL)
		logError() << "Could not open fault checkpoint file";

	// Read and broadcast header
	checkMPIErr(setHeaderView(file));

	Header header;
	if (rank() == 0)
		checkMPIErr(MPI_File_read(file, &header, 1, headerType(), MPI_STATUS_IGNORE));

	MPI_Bcast(&header, 1, headerType(), 0, comm());
	timestepFault = header.timestepFault;

	double* data[NUM_VARIABLES] = {mu, slipRate1, slipRate2, slip, slip1, slip2, state, strength};

	// Read data
	checkMPIErr(setDataView(file));
	for (unsigned int i = 0; i < NUM_VARIABLES; i++)
		checkMPIErr(MPI_File_read_all(file, data[i], numSides() * numBndGP(), MPI_DOUBLE, MPI_STATUS_IGNORE));

	// Close the file
	checkMPIErr(MPI_File_close(&file));
}

void seissol::checkpoint::mpio::Fault::write(int timestepFault)
{
	EPIK_TRACER("CheckPointFault_write");
	SCOREP_USER_REGION("CheckPointFault_write", SCOREP_USER_REGION_TYPE_FUNCTION);

	if (numSides() == 0)
		return;

	logInfo(rank()) << "Checkpoint backend: Writing fault.";

	// Write the header
	writeHeader(timestepFault);

	// Save data
	EPIK_USER_REG(r_write_wavefield, "checkpoint_write_fault");
	SCOREP_USER_REGION_DEFINE(r_write_fault);
	EPIK_USER_START(r_write_wavefield);
	SCOREP_USER_REGION_BEGIN(r_write_fault, "checkpoint_write_fault", SCOREP_USER_REGION_TYPE_COMMON);

	checkMPIErr(setDataView(file()));

	for (unsigned int i = 0; i < NUM_VARIABLES; i++)
		checkMPIErr(MPI_File_write_all(file(), const_cast<double*>(data(i)), numSides() * numBndGP(), MPI_DOUBLE, MPI_STATUS_IGNORE));

	EPIK_USER_END(r_write_fault);
	SCOREP_USER_REGION_END(r_write_fault);

	// Finalize the checkpoint
	finalizeCheckpoint();

	logInfo(rank()) << "Checkpoint backend: Writing fault. Done.";
}

bool seissol::checkpoint::mpio::Fault::validate(MPI_File file)
{
	int result = true;

	if (rank() == 0) {
		Header header;

		// Check the header
		MPI_File_read(file, &header, 1, headerType(), MPI_STATUS_IGNORE);

		if (header.identifier != identifier()) {
			logWarning() << "Checkpoint identifier does match";
			result = false;
		}
	}

	// Make sure everybody knows the result of the validation
	MPI_Bcast(&result, 1, MPI_INT, 0, comm());

	return result;
}

void seissol::checkpoint::mpio::Fault::writeHeader(int timestepFault)
{
	EPIK_TRACER("checkpoint_write_fault_header");
	SCOREP_USER_REGION("checkpoint_write_fault_header", SCOREP_USER_REGION_TYPE_FUNCTION);

	checkMPIErr(setHeaderView(file()));

	if (rank() == 0) {
		Header header;
		header.identifier = identifier();
		header.timestepFault = timestepFault;

		checkMPIErr(MPI_File_write(file(), &header, 1, headerType(), MPI_STATUS_IGNORE));
	}
}
