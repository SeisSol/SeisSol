/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2014-2015, SeisSol Group
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
 * Main C++ SeisSol file
 */

#include "SeisSol.h"

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "utils/timeutils.h"

void seissol::SeisSol::init(int argc, char* argv[])
{
	m_mpi.init(argc, argv);

	const int rank = m_mpi.rank();

  // Print welcome message
  logInfo(rank) << "Welcome to SeisSol";
  logInfo(rank) << "Copyright (c) 2012-2015, SeisSol Group";
  logInfo(rank) << "Built on:" << __DATE__ << __TIME__
#ifdef GENERATEDKERNELS
                << "(generated kernels)"
#endif // GENERATEDKERNELS
  ;
#ifdef _OPENMP
  logInfo(rank) << "Version:" << SEISSOL_VERSION_STRING;
  logInfo(rank) << "Using OMP with #threads/rank:" << omp_get_max_threads();
#ifdef USE_MPI
#ifdef USE_COMM_THREAD
  logWarning(rank) << "Running with communication thread in hybrid mode!";
  logWarning(rank) << "Make sure that OMP is just using N-1 cores!";
  logWarning(rank) << "The communication thread will pin itself to OS core-id N-1!";
  logWarning(rank) << "Check the pinning settings of your OMP implementation!";
  logWarning(rank) << "Do not pin OMP threads to HW threads of the last OS core!";
  logWarning(rank) << "Do not use more than 1 MPI rank per node!";
#endif
#endif
#endif // _OPENMP
}

void seissol::SeisSol::finalize()
{
	const int rank = m_mpi.rank();

	m_mpi.finalize();

	logInfo(rank) << "SeisSol done. Goodbye.";
}

seissol::SeisSol seissol::SeisSol::main;

// Fortran main definition
extern "C" {

void fortran_main();

}

// Start here

int main(int argc, char* argv[])
{
	EPIK_TRACER("SeisSol");
	SCOREP_USER_REGION("SeisSol", SCOREP_USER_REGION_TYPE_FUNCTION);

	// Initialize SeisSol
	seissol::SeisSol::main.init(argc, argv);

	// Initialize Fortan Part and run SeisSol
	fortran_main();

	// Finalize SeisSol
	seissol::SeisSol::main.finalize();
}
