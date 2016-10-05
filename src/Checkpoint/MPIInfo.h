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

#ifndef CHECKPOINT_MPIINFO_H
#define CHECKPOINT_MPIINFO_H

#include <mpi.h>

#include <string>

#include <utils/env.h>
#include <utils/stringutils.h>

namespace seissol
{

namespace checkpoint
{

/**
 * Generates the MPI_INFO object for MPI-IO based checkpoints
 */
class MPIInfo
{
private:
	/** The MPI_INFO object */
	MPI_Info m_info;

public:
	MPIInfo()
	{
		MPI_Info_create(&m_info);

		// Set parameters
		const int numParams = 5;
		const char* params[numParams] = {"cb_nodes",
			"romio_cb_read", "romio_cb_write",
			"romio_ds_read", "romio_ds_write"};

		for (int i = 0; i < numParams; i++) {
			std::string param(params[i]);
			utils::StringUtils::toUpper(param);
			std::string envName = "SEISSOL_CHECKPOINT_" + param;

			const char* value = utils::Env::get<const char*>(envName.c_str(), 0L);
			if (value) {
				MPI_Info_set(m_info, const_cast<char*>(params[i]), const_cast<char*>(value));
			}
		}
	}

	~MPIInfo()
	{
		MPI_Info_free(&m_info);
	}

	MPI_Info get() const
	{
		return m_info;
	}
};

}

}

#endif // CHECKPOINT_MPIINFO_H
