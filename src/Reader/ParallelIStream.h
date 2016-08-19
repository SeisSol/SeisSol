/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2014-2016, SeisSol Group
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
 */

#ifndef PARALLEL_ISTREAM_H
#define PARALLEL_ISTREAM_H

#include "Parallel/MPI.h"

#include "utils/logger.h"

#include <fstream>
#include <sstream>
#include <vector>

class ParallelIStream : public std::wistringstream
{
private:
	/** Contains the content of the file */
	std::vector<wchar_t> m_buffer;

public:
	ParallelIStream()
	{
	}

	ParallelIStream(const char* filename)
	{
		open(filename);
	}

	void open(const char* filename)
	{
#ifdef USE_MPI
		MPI_Comm comm = seissol::MPI::mpi.comm();
		const int rank = seissol::MPI::mpi.rank();

		if (rank == 0) {
#endif // USE_MPI

			// Open the file and set the buffer
			// This only done on rank 0 when running in parallel
			std::wifstream file(filename);
			if (!file)
				logError() << "Could not open file" << filename;

			file.seekg(0, std::ios::end);
			m_buffer.resize(file.tellg());
			file.seekg(0, std::ios::beg);

			file.read(&m_buffer[0], m_buffer.size());

			rdbuf()->pubsetbuf(&m_buffer[0], m_buffer.size());

#ifdef USE_MPI
			// Broadcast the size and the content of the file
			unsigned long size = m_buffer.size();
			MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG, 0, comm);
			MPI_Bcast(&m_buffer[0], size, MPI_WCHAR, 0, comm);
		} else {
			unsigned long size;
			MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG, 0, comm);

			m_buffer.resize(size);
			MPI_Bcast(&m_buffer[0], size, MPI_WCHAR, 0, comm);

			rdbuf()->pubsetbuf(&m_buffer[0], size);
		}
#endif // USE_MPI
	}
};

#endif // PARALLEL_ISTREAM_H
