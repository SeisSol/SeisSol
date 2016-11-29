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
 * Velocity field reader Fortran interface
 */

#include "Parallel/MPI.h"

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "utils/logger.h"

#include "AsagiReader.h"
#include "Modules/Modules.h"
#include "Monitoring/instrumentation.fpp"
#include "Monitoring/Stopwatch.h"

namespace seissol
{

struct Setter
{
	/** The material values we need to write to */
	double* const materialValues;

	/** Number of cells */
	const unsigned int numElements;

	Setter(double* materialValues, unsigned int numElements)
		: materialValues(materialValues), numElements(numElements)
	{
	}


	void set(unsigned int i, const float* material)
	{
		materialValues[i] = material[0];
		materialValues[i+numElements] = material[1];
		materialValues[i+2*numElements] = material[2];
	}
};

}

extern "C"
{

/**
 * Reads a velocity field using ASAGI
 *
 * @warning Since the ASAGI code might use threads as well as OpenMP, we need to switch
 *  instrumentation off for some parts. This can only be done with Scalasca 2.x.
 *  Therefore, ASAGI cannot run with Scalasca 1.x.
 */
void read_velocity_field(const char* file, int numElements, const seissol::asagi::vertex_t* baryCenters,
        double defaultRho, double defaultMu, double defaultLambda, double* materialValues)
{
#ifdef USE_ASAGI
	SCOREP_USER_REGION("read_velocity_field", SCOREP_USER_REGION_TYPE_FUNCTION);

	const int rank = seissol::MPI::mpi.rank();

	seissol::asagi::AsagiReader reader("SEISSOL_ASAGI");

	// Use manual timing since Score-P does not work
	Stopwatch stopwatch;
	stopwatch.start();

	logInfo(rank) << "Initializing velocity field.";
	reader.open(file, 3);

	double time = stopwatch.stop();
	logInfo(rank) << "Velocity field opened in" << time << "sec.";

	// Start reading the data
	stopwatch.start();

	seissol::Setter setter(materialValues, numElements);
	float defaultValues[3] = {defaultRho, defaultMu, defaultLambda};
	unsigned long outside;
	reader.read(baryCenters, numElements,
		setter, defaultValues,
		outside);
	if (outside > 0)
		logWarning(rank) << "Found" << outside << "cells of the given velocity field.";

	time = stopwatch.stop();
	logInfo(rank) << "Velocity field initialized in" << time << "sec.";

	logInfo(rank) << "Initializing velocity field. Done.";

#else // USE_ASAGI
	logError() << "This version does not support ASAGI";
#endif // USE_ASAGI
}

}
