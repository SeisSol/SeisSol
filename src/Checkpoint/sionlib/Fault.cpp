/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Gilbert Brietzke (gilbert.brietzke AT lrz.de, http://www.lrz.de)
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

#include "Fault.h"
#include "Kernels/precision.hpp"

bool seissol::checkpoint::sionlib::Fault::init(unsigned int numSides, unsigned int numBndGP,
		unsigned int groupSize)
{
	if (groupSize != 1)
		// TODO To read the sionlib file, we must use the same number of processes
		// (otherwise it gets very complicated).
		logError() << "The SIONlib backend does not support asynchronous MPI mode yet.";

	setChunkElementCount(numSides * numBndGP);

	seissol::checkpoint::Fault::init(numSides, numBndGP, groupSize);

	if (numSides == 0)
		return true;

	return exists();
}

void seissol::checkpoint::sionlib::Fault::load(int &timestepFault, double* mu, double* slipRate1, double* slipRate2,
	double* slip, double* slip1, double* slip2, double* state, double* strength)
{
	if (numSides() == 0)
		return;

	logInfo(rank()) << "Loading fault checkpoint";

	seissol::checkpoint::CheckPoint::setLoaded();

	int file = open(linkFile(), readMode());
	checkErr(file);

	// Read identifier
	unsigned long id;
	checkErr(sion_coll_fread(&id, sizeof(id), 1, file), 1);

	// Read header
	checkErr(sion_coll_fread(&timestepFault, sizeof(timestepFault), 1, file), 1);

	double* data[NUM_VARIABLES] = {mu, slipRate1, slipRate2, slip, slip1, slip2, state, strength};

	// Read data
	for (unsigned int i = 0; i < NUM_VARIABLES; i++)
		checkErr(sion_coll_fread(data[i], sizeof(real), numSides()*numBndGP(), file),
				numSides()*numBndGP());

	// Close the file
	sionClose(file);
}

void seissol::checkpoint::sionlib::Fault::write(int timestepFault)
{
	SCOREP_USER_REGION("CheckPointFault_write", SCOREP_USER_REGION_TYPE_FUNCTION);

	if (numSides() == 0)
		return;

	logInfo(rank()) << "Checkpoint backend: Writing fault.";

	int file = open(dataFile(odd()), writeMode());
	checkErr(file);

	// Write the header
	SCOREP_USER_REGION_DEFINE(r_write_header);
	SCOREP_USER_REGION_BEGIN(r_write_header, "checkpoint_write_fault_header", SCOREP_USER_REGION_TYPE_COMMON);

	unsigned long id = identifier();
	checkErr(sion_coll_fwrite(&id, sizeof(id), 1, file), 1);
	checkErr(sion_coll_fwrite(&timestepFault, sizeof(timestepFault), 1, file), 1);

	SCOREP_USER_REGION_END(r_write_header);

	// Save data
	SCOREP_USER_REGION_DEFINE(r_write_fault);
	SCOREP_USER_REGION_BEGIN(r_write_fault, "checkpoint_write_fault", SCOREP_USER_REGION_TYPE_COMMON);

	for (unsigned int i = 0; i < NUM_VARIABLES; i++)
	    checkErr(sion_coll_fwrite(data(i), sizeof(real), numSides() * numBndGP(), file),
	    		numSides() * numBndGP());

	SCOREP_USER_REGION_END(r_write_fault);

	// Finalize the checkpoint
	finalizeCheckpoint(file);

	logInfo(rank()) << "Checkpoint backend: Writing fault. Done.";
}
