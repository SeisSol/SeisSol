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

#include "Wavefield.h"

bool seissol::checkpoint::posix::Wavefield::init(size_t headerSize, unsigned long numDofs, unsigned int groupSize)
{
	seissol::checkpoint::Wavefield::init(headerSize, numDofs, groupSize);

	return exists();
}

void seissol::checkpoint::posix::Wavefield::load(real* dofs)
{
	logInfo(rank()) << "Loading wave field checkpoint";

	seissol::checkpoint::CheckPoint::setLoaded();

	int file = open();
	checkErr(file);

	// Read header
	checkErr(read(file, header().data(), header().size()), header().size());

	// Skip other processes before this in the group
	checkErr(lseek64(file, groupOffset() * sizeof(real), SEEK_CUR));

	// Convert to char* to do pointer arithmetic
	char* buffer = reinterpret_cast<char*>(dofs);
	unsigned long left = numDofs()*sizeof(real);

	// Read dofs
	while (left > 0) {
		unsigned long readSize = read(file, buffer, left);
		if (readSize <= 0)
			checkErr(readSize, left);
		buffer += readSize;
		left -= readSize;
	}

	// Close the file
	checkErr(::close(file));
}

void seissol::checkpoint::posix::Wavefield::write(const void* header, size_t headerSize)
{
	EPIK_TRACER("CheckPoint_write");
	SCOREP_USER_REGION("CheckPoint_write", SCOREP_USER_REGION_TYPE_FUNCTION);

	logInfo(rank()) << "Checkpoint backend: Writing.";

	// Start at the beginning
	checkErr(lseek64(file(), 0, SEEK_SET));

	// Write the header
	EPIK_USER_REG(r_write_header, "checkpoint_write_header");
	SCOREP_USER_REGION_DEFINE(r_write_header);
	EPIK_USER_START(r_write_header);
	SCOREP_USER_REGION_BEGIN(r_write_header, "checkpoint_write_header", SCOREP_USER_REGION_TYPE_COMMON);

	checkErr(::write(file(), header, headerSize), headerSize);

	EPIK_USER_END(r_write_header);
	SCOREP_USER_REGION_END(r_write_header);

	// Save data
	EPIK_USER_REG(r_write_wavefield, "checkpoint_write_wavefield");
	SCOREP_USER_REGION_DEFINE(r_write_wavefield);
	EPIK_USER_START(r_write_wavefield);
	SCOREP_USER_REGION_BEGIN(r_write_wavefield, "checkpoint_write_wavefield", SCOREP_USER_REGION_TYPE_COMMON);

	// Convert to char* to do pointer arithmetic
	const char* buffer = reinterpret_cast<const char*>(dofs());
	unsigned long left = numDofs()*sizeof(real);
	if (alignment()) {
		left = (left + alignment() - 1) / alignment();
		left *= alignment();
	}

	while (left > 0) {
		unsigned long written = ::write(file(), buffer, left);
		if (written <= 0)
			checkErr(written, left);
		buffer += written;
		left -= written;
	}

	EPIK_USER_END(r_write_wavefield);
	SCOREP_USER_REGION_END(r_write_wavefield);

	// Finalize the checkpoint
	finalizeCheckpoint();

	logInfo(rank()) << "Checkpoint backend: Writing. Done.";
}
