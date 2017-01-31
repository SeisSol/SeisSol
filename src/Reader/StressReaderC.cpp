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

#include "Parallel/MPI.h"

#include <cstring>

#include "utils/env.h"
#include "utils/logger.h"

#include "AsagiReader.h"
#include "Monitoring/Stopwatch.h"

#ifdef USE_ASAGI
/** The ASAGI reader */
static seissol::asagi::AsagiReader reader("SEISSOL_ASAGI_STRESS");

/** The number of variables found in ASAGI */
unsigned int numVariables;

/** Counter for values outside the box */
static unsigned long outside = 0;

/** The stopwatch for timing */
static Stopwatch stopwatch;
#endif // USE_ASAGI

namespace seissol
{

namespace stress_reader
{

struct Setter
{
	const unsigned int numValues;

	/** Destination buffer for stress values */
	double * const buffer;

	Setter(unsigned int numValues, double * const buffer)
		: numValues(numValues), buffer(buffer)
	{
	}

	void set(const float* stressValues)
	{
		for (unsigned int i = 0; i < 6; i++)
			buffer[i] = stressValues[i];
		if (numValues > 6) {
			// Substract the p from sxx, syy, szz
			buffer[0] -= stressValues[6];
			buffer[1] -= stressValues[6];
			buffer[2] -= stressValues[6];
		}
	}
};

}

}

extern "C"
{

/**
 * @warning Instrumentation does not work with ASAGI. We use the stopwatch to messure time.
 *
 * @todo Support stress initialization for the whole domain
 */
void open_stress_field(const char* file)
{
#ifdef USE_ASAGI
	const int rank = seissol::MPI::mpi.fault.rank();

	// Use manual timing since Score-P does not work
	stopwatch.start();

	bool sparse = utils::Env::get<bool>("SEISSOL_ASAGI_STRESS_SPARSE", false);

	logInfo(rank) << "Initializing stress field.";
	numVariables = reader.open(file, sparse
#ifdef USE_MPI
		, seissol::MPI::mpi.fault.comm()
#endif // USE_MPI
	);
	switch (numVariables) {
	case 6:
	case 7:
		break;
	default:
		logError() << "Invalid number of variables in stress input";
	}

	double time = stopwatch.stop();
	logInfo(rank) << "Stress field opened in" << time << "sec.";

	// Start reading the data
	stopwatch.start();

#else // USE_ASAGI
	logError() << "This version does not support ASAGI";
#endif // USE_ASAGI
}

/**
 * @warning Instrumentation does not work with ASAGI. We use the stopwatch to messure time.
 */
void read_stress(double x, double y, double z, double values[6], const float defaultValues[7])
{
#ifdef USE_ASAGI
	const seissol::asagi::vertex_t coord = {x, y, z};

	seissol::stress_reader::Setter setter(numVariables, values);
	outside += reader.readValue(coord, setter, defaultValues);

#else // USE_ASAGI
	logError() << "This version does not support ASAGI";
#endif // USE_ASAGI
}

/**
 * @warning Instrumentation does not work with ASAGI. We use the stopwatch to messure time.
 *
 * @todo Support stress initialization for the whole domain
 */
void close_stress_field()
{
#ifdef USE_ASAGI
	const int rank = seissol::MPI::mpi.fault.rank();

#ifdef USE_MPI
	if (rank == 0)
		MPI_Reduce(MPI_IN_PLACE, &outside, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, seissol::MPI::mpi.fault.comm());
	else
		MPI_Reduce(&outside, 0L, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, seissol::MPI::mpi.fault.comm());
#endif // USE_MPI
	if (outside > 0)
		logWarning(rank) << "Found" << outside << "cells outside of the given stress field.";

	double time = stopwatch.stop();
	logInfo(rank) << "Stress field initialized in" << time << "sec.";

	reader.close();

	logInfo(rank) << "Initializing stress field. Done.";
#else // USE_ASAGI
	logError() << "This version does not support ASAGI";
#endif // USE_ASAGI
}

}
