/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Asynchronous I/O
 */

#ifndef ASYNCIO_H
#define ASYNCIO_H

#include "Parallel/MPI.h"

#include <algorithm>

#include "utils/env.h"
#include "utils/logger.h"

#include "async/Dispatcher.h"

namespace seissol
{

namespace io
{

class AsyncIO : public async::Dispatcher
{
public:
	/**
	 * @return False if this rank is an MPI executor that does not contribute to the
	 *  computation.
	 */
	bool init()
	{
		async::Dispatcher::init();

#ifdef USE_MPI
		seissol::MPI::mpi.setComm(commWorld());
		// TODO Update fault communicator (not really sure how we can do this at this point)
#endif // USE_MPI

		return dispatch();
	}

	void finalize()
	{
		// Call parent class
		async::Dispatcher::finalize();

#ifdef USE_MPI
		// Reset the MPI communicator
		seissol::MPI::mpi.setComm(MPI_COMM_WORLD);
#endif // USE_MPI
	}
};

}

}

#endif // ASYNCIO_H
