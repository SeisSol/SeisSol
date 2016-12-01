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
 * MPI Wrapper
 */

#ifndef MPI_H
#define MPI_H

#ifndef USE_MPI
#include "MPIDummy.h"
#else // USE_MPI

#include <mpi.h>

#include "utils/logger.h"

#include "MPIBasic.h"

#endif // USE_MPI

namespace seissol
{

#ifndef USE_MPI
typedef MPIDummy MPI;
#else // USE_MPI

/**
 * MPI handling.
 *
 * Make sure only one instance of this class exists!
 */
class MPI : public MPIBasic
{
private:
	MPI_Comm m_comm;

private:
	MPI()
		: m_comm(MPI_COMM_NULL)
	{ }

public:
	~MPI()
	{ }

	/**
	 * Initialize MPI
	 */
	void init(int &argc, char** &argv)
	{
		int required = (m_threadsafe ? MPI_THREAD_MULTIPLE : MPI_THREAD_SINGLE);
		int provided;
		MPI_Init_thread(&argc, &argv, required, &provided);

		setComm(MPI_COMM_WORLD);

		// Test this after setComm() to get the correct m_rank
		if (required < provided)
			logWarning(m_rank) << utils::nospace << "Required MPI thread support (" << required
				<< ") is smaller than provided thread support (" << provided << ").";
	}

	void setComm(MPI_Comm comm)
	{
		m_comm = comm;

		MPI_Comm_rank(comm, &m_rank);
		MPI_Comm_size(comm, &m_size);
	}

	/**
	 * @return The main communicator for the application
	 */
	MPI_Comm comm() const
	{
		return m_comm;
	}

	void barrier(MPI_Comm comm) const
	{
		MPI_Barrier(comm);
	}

	/**
	 * Finalize MPI
	 */
	void finalize()
	{
		fault.finalize();

		MPI_Finalize();
	}

public:
	/** The only instance of the class */
	static MPI mpi;
};

#endif // USE_MPI

}

#endif // MPI_H
