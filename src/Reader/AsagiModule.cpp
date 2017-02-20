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
 * Velocity field reader Fortran interface
 */

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include <string>

#include "utils/env.h"

#include "AsagiModule.h"

seissol::asagi::AsagiModule::AsagiModule()
	: m_mpiMode(getMPIMode()), m_totalThreads(getTotalThreads())
{
	// Register for the pre MPI hook
	Modules::registerHook(*this, seissol::PRE_MPI);

	// Emit a warning/error later
	// TODO use a general logger that can buffer log messages and emit them later
	if (m_mpiMode == MPI_UNKNOWN) {
		Modules::registerHook(*this, seissol::POST_MPI_INIT);
	} else if (m_mpiMode == MPI_COMM_THREAD && m_totalThreads == 1) {
		m_mpiMode = MPI_WINDOWS;

		Modules::registerHook(*this, seissol::POST_MPI_INIT);
	}
}

seissol::asagi::AsagiModule seissol::asagi::AsagiModule::instance;

seissol::asagi::MPI_Mode seissol::asagi::AsagiModule::getMPIMode()
{
#ifdef USE_MPI
	std::string mpiModeName = utils::Env::get(ENV_MPI_MODE, "WINDOWS");
	if (mpiModeName == "WINDOWS")
		return MPI_WINDOWS;
	if (mpiModeName == "COMM_THREAD")
		return MPI_COMM_THREAD;
	if (mpiModeName == "OFF")
		return MPI_OFF;

	return MPI_UNKNOWN;
#else // USE_MPI
	return MPI_OFF;
#endif // USE_MPI
}

int seissol::asagi::AsagiModule::getTotalThreads()
{
	int totalThreads = 1;

#ifdef _OPENMP
	totalThreads = omp_get_max_threads();
#ifdef USE_COMM_THREAD
	totalThreads++;
#endif // USE_COMM_THREAD
#endif // _OPENMP

	return totalThreads;
}

const char* seissol::asagi::AsagiModule::ENV_MPI_MODE = "SEISSOL_ASAGI_MPI_MODE";