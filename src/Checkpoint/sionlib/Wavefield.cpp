/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Gilbert Brietzke (gilbert.brietzke AT LRZ.de, http://www.lrz.de)
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

#include "Wavefield.h"

bool seissol::checkpoint::sionlib::Wavefield::init(size_t headerSize, unsigned long numDofs, unsigned int groupSize)
{
	if (groupSize != 1)
		// TODO To read the sionlib file, we must use the same number of processes
		// (otherwise it gets very complicated).
		logError() << "The SIONlib backend does not support asynchronous MPI mode yet.";

	setChunkElementCount(numDofs);

	seissol::checkpoint::Wavefield::init(headerSize, numDofs, groupSize);

	logInfo(rank()) << "Using SIONlib mode" << writeMode() << "for writing checkpoints";

	return exists();
}

void seissol::checkpoint::sionlib::Wavefield::load(real* dofs)
{
	logInfo(rank()) << "Loading wave field checkpoint";

	seissol::checkpoint::CheckPoint::setLoaded();

	int file = open(linkFile(), readMode());
	checkErr(file);

	// Read header
	checkErr(sion_coll_fread(header().data(), header().size(), 1, file), 1);

	// Read dofs
	checkErr(sion_coll_fread(dofs, sizeof(double), numDofs(), file), numDofs());

	// Close the file
	sionClose(file);
}

void seissol::checkpoint::sionlib::Wavefield::write(const void* header, size_t headerSize)
{
	SCOREP_USER_REGION("CheckPoint_write", SCOREP_USER_REGION_TYPE_FUNCTION);

	logInfo(rank()) << "Checkpoint backend: Writing.";

	int file = open(dataFile(odd()), writeMode());
	checkErr(file);

	// Write the header
	SCOREP_USER_REGION_DEFINE(r_write_header);
	SCOREP_USER_REGION_BEGIN(r_write_header, "checkpoint_write_header", SCOREP_USER_REGION_TYPE_COMMON);

	checkErr(sion_coll_fwrite(header, headerSize, 1, file), 1);

	SCOREP_USER_REGION_END(r_write_header);

	// Save data
	SCOREP_USER_REGION_DEFINE(r_write_wavefield);
	SCOREP_USER_REGION_BEGIN(r_write_wavefield, "checkpoint_write_wavefield", SCOREP_USER_REGION_TYPE_COMMON);

	checkErr(sion_coll_fwrite(dofs(), sizeof(real), numDofs(), file), numDofs());

	SCOREP_USER_REGION_END(r_write_wavefield);

	// Finalize the checkpoint
	finalizeCheckpoint(file);

	logInfo(rank()) << "Checkpoint backend: Writing. Done.";
}

