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

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include <string>

#ifdef USE_ASAGI
#include <asagi.h>
#endif // USE_ASAGI

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

#include "utils/env.h"
#include "utils/logger.h"

#include "glm/vec3.hpp"
#include "glm/vector_relational.hpp"
#include "glm/gtc/round.hpp"

#include "Parallel/MPI.h"
#include "Monitoring/instrumentation.fpp"
#include "Monitoring/Stopwatch.h"

typedef double vertex_t[3];

enum MPI_Mode
{
	MPI_OFF, MPI_WINDOWS, MPI_COMM_THREAD
};

enum NUMACache_Mode
{
	NUMA_OFF, NUMA_ON, NUMA_CACHE
};

#ifdef USE_ASAGI
static int getTotalThreads()
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

/**
 * Get the MPI mode from the environment
 */
static MPI_Mode getMPIMode(int totalThreads)
{
	const int rank = seissol::MPI::mpi.rank();

#ifdef USE_MPI
	const char* mpiModeName = utils::Env::get("SEISSOL_ASAGI_MPI_MODE", "WINDOWS");
	if (strcmp(mpiModeName, "WINDOWS") == 0)
		return MPI_WINDOWS;
	if (strcmp(mpiModeName, "COMM_THREAD") == 0) {
		if (totalThreads > 1)
			return MPI_COMM_THREAD;

		logWarning(rank) << "Running with only one OMP thread."
				<< "Using MPI window communication instead of threads.";
		return MPI_WINDOWS;
	}
	if (strcmp(mpiModeName, "OFF") == 0)
		return MPI_OFF;

	logError() << "Unknown MPI mode:" << mpiModeName;
#endif // USE_MPI
	return MPI_OFF;
}

static NUMACache_Mode getNUMAMode()
{
	const char* numaModeName = utils::Env::get("SEISSOL_ASAGI_NUMA_MODE", "ON");

	if (strcmp(numaModeName, "ON") == 0)
		return NUMA_ON;
	if (strcmp(numaModeName, "OFF") == 0)
		return NUMA_OFF;
	if (strcmp(numaModeName, "CACHE") == 0)
		return NUMA_CACHE;

	logError() << "Unknown NUMA mode:" << numaModeName;
	return NUMA_OFF;
}
#endif // USE_ASAGI

extern "C"
{

/**
 * Reads a velocity field using ASAGI
 *
 * @warning Since the ASAGI code might use threads as well as OpenMP, we need to switch
 *  instrumentation off for some parts. This can only be done with Scalasca 2.x.
 *  Therefore, ASAGI cannot run with Scalasca 1.x.
 */
void read_velocity_field(const char* file, int numElements, const vertex_t* baryCenters,
        double defaultRho, double defaultMu, double defaultLambda, double* materialValues)
{
#ifdef USE_ASAGI
	SCOREP_USER_REGION("read_velocity_field", SCOREP_USER_REGION_TYPE_FUNCTION);

	// Use manual timing since Score-P does not work
	Stopwatch stopwatch;
	stopwatch.start();

	const int rank = seissol::MPI::mpi.rank();

	logInfo(rank) << "Initializing velocity field.";

	SCOREP_USER_REGION_DEFINE(r_asagi_init);
	SCOREP_USER_REGION_BEGIN(r_asagi_init, "asagi_init", SCOREP_USER_REGION_TYPE_COMMON);

	asagi::Grid* grid = asagi::Grid::createArray();

	int totalThreads = getTotalThreads();

	// Set MPI mode
	MPI_Mode mpiMode = getMPIMode(totalThreads);
#ifdef USE_MPI
	if (mpiMode != MPI_OFF) {
		grid->setComm(seissol::MPI::mpi.comm());

		if (mpiMode == MPI_COMM_THREAD) {
			// WARNING: This assumes that if we use the communication thread,
			// only one process per rank is started!
			//SCOREP_RECORDING_OFF();
			asagi::Grid::startCommThread();
			//SCOREP_RECORDING_ON();

			grid->setParam("MPI_COMMUNICATION", "THREAD");
		}
	}
#endif // USE_MPI

	// Set NUMA mode
	unsigned int asagiThreads = utils::Env::get("SEISSOL_ASAGI_NUM_THREADS", 0u);
	if (asagiThreads == 0)
		asagiThreads = totalThreads;
	else if (asagiThreads > totalThreads) {
		logWarning(rank) << "Only" << totalThreads
				<< "threads can be used for ASAGI initialization.";
		asagiThreads = totalThreads;
	}

	if (mpiMode == MPI_COMM_THREAD)
		asagiThreads--; // one thread is used for communication

	grid->setThreads(asagiThreads);

	switch (getNUMAMode()) {
	case NUMA_ON:
		grid->setParam("NUMA_COMMUNICATION", "ON");
		break;
	case NUMA_OFF:
		grid->setParam("NUMA_COMMUNICATION", "OFF");
		break;
	case NUMA_CACHE:
		grid->setParam("NUMA_COMMUNICATION", "CACHE");
		break;
	}

	// Set vertex centered grid
	grid->setParam("VALUE_POSITION", "VERTEX_CENTERED");

	// Set additional parameters
	std::string blockSize = utils::Env::get("SEISSOL_ASAGI_BLOCK_SIZE", "64");
	grid->setParam("BLOCK_SIZE_0", blockSize.c_str());
	grid->setParam("BLOCK_SIZE_1", blockSize.c_str());
	grid->setParam("BLOCK_SIZE_2", blockSize.c_str());

	std::string cacheSize = utils::Env::get("SEISSOL_ASAGI_CACHE_SIZE", "128");
	grid->setParam("CACHE_SIZE", cacheSize.c_str());

	grid->setParam("VARIABLE", "data");

	// Read the data
	//SCOREP_RECORDING_OFF();
#ifdef _OPENMP
	#pragma omp parallel num_threads(asagiThreads)
#endif // _OPENMP
	{
		asagi::Grid::Error err = grid->open(file);
		if (err != asagi::Grid::SUCCESS)
			logError() << "Could not open ASAGI grid:" << err;
	}
	//SCOREP_RECORDING_ON();

	if (grid->getVarSize() != 3*sizeof(float))
		logError() << "Invalid variable size in material file";

	double time = stopwatch.stop();
	logInfo(rank) << "Velocity field opened in" << time << "sec.";

	SCOREP_USER_REGION_END(r_asagi_init);

	// Grid dimensions
	const glm::dvec3 min(grid->getMin(0), grid->getMin(1), grid->getMin(2));
	const glm::dvec3 max(grid->getMax(0), grid->getMax(1), grid->getMax(2));
	const glm::dvec3 delta(grid->getDelta(0), grid->getDelta(1), grid->getDelta(2));

	SCOREP_USER_REGION_DEFINE(r_asagi_read);
	SCOREP_USER_REGION_BEGIN(r_asagi_read, "asagi_read", SCOREP_USER_REGION_TYPE_COMMON);

	stopwatch.start();

	// Initialize the values in SeisSol
	unsigned long outside = 0;

	//SCOREP_RECORDING_OFF();
#ifdef _OPENMP
	#pragma omp parallel for num_threads(asagiThreads) reduction(+: outside)
#endif // _OPENMP
	for (int i = 0; i < numElements; i++) {
		const glm::dvec3 baryCenter(baryCenters[i][0], baryCenters[i][1], baryCenters[i][2]);

		// Compute the surrounding coordinates
		glm::dvec3 shiftedBary = baryCenter - min;
		glm::dvec3 lowCoord = glm::floorMultiple(shiftedBary, delta) + min;
		glm::dvec3 highCoord = glm::ceilMultiple(shiftedBary, delta) + min;

		// Fix low/high if they are outside of the domain (but the bary center is inside)
		// -> Should not be necessary if we use vertex centered grid
//		if (glm::all(glm::greaterThanEqual(baryCenter, min))
//			&& glm::any(glm::lessThan(lowCoord, min)))
//			lowCoord = glm::max(lowCoord, min);
//		if (glm::all(glm::lessThanEqual(baryCenter, max))
//			&& glm::any(glm::greaterThan(highCoord, max)))
//			highCoord = glm::min(highCoord, max);

		if (glm::any(glm::lessThan(lowCoord, min))
			|| glm::any(glm::greaterThan(highCoord, max))) {
			// Outside the ASAGI domain

			materialValues[i] = defaultRho;
			materialValues[i+numElements] = defaultMu;
			materialValues[i+2*numElements] = defaultLambda;

			outside++;
		} else {
			// Get all 8 material values
			float material[8][3];
			double coord[3];
			coord[0] = lowCoord.x; coord[1] = lowCoord.y; coord[2] = lowCoord.z;
			grid->getBuf(material[0], coord);
			coord[0] = lowCoord.x; coord[1] = lowCoord.y; coord[2] = highCoord.z;
			grid->getBuf(material[1], coord);
			coord[0] = lowCoord.x; coord[1] = highCoord.y; coord[2] = lowCoord.z;
			grid->getBuf(material[2], coord);
			coord[0] = lowCoord.x; coord[1] = highCoord.y; coord[2] = highCoord.z;
			grid->getBuf(material[3], coord);
			coord[0] = highCoord.x; coord[1] = lowCoord.y; coord[2] = lowCoord.z;
			grid->getBuf(material[4], coord);
			coord[0] = highCoord.x; coord[1] = lowCoord.y; coord[2] = highCoord.z;
			grid->getBuf(material[5], coord);
			coord[0] = highCoord.x; coord[1] = highCoord.y; coord[2] = lowCoord.z;
			grid->getBuf(material[6], coord);
			coord[0] = highCoord.x; coord[1] = highCoord.y; coord[2] = highCoord.z;
			grid->getBuf(material[7], coord);

			// Do trilinear interpolation:
			// https://en.wikipedia.org/wiki/Trilinear_interpolation
			glm::dvec3 d = (baryCenter - lowCoord) / delta;

			double interpolMaterial[3];
			for (unsigned int j = 0; j < 3; j++) {
				double c0[4];
				for (unsigned int k = 0; k < 4; k++)
					c0[k] = material[k][j]*(1-d.x) + material[k+4][j]*d.x;

				double c1[2];
				for (unsigned int k = 0; k < 2; k++)
					c1[k] = c0[k]*(1-d.y) + c0[k+1]*d.y;

				interpolMaterial[j] = c1[0]*(1-d.z) + c1[1]*d.z;
			}

			materialValues[i] = interpolMaterial[0];
			materialValues[i+numElements] = interpolMaterial[1];
			materialValues[i+2*numElements] = interpolMaterial[2];
		}
	}
	//SCOREP_RECORDING_ON();

#ifdef USE_MPI
	if (rank == 0)
		MPI_Reduce(MPI_IN_PLACE, &outside, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, seissol::MPI::mpi.comm());
	else
		MPI_Reduce(&outside, 0L, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, seissol::MPI::mpi.comm());
#endif // USE_MPI
	if (outside > 0)
		logWarning(rank) << "Found" << outside << "cells of the given velocity field.";

	time = stopwatch.stop();
	logInfo(rank) << "Velocity field initialized in" << time << "sec.";

	SCOREP_USER_REGION_END(r_asagi_read);

	// Cleanup
	delete grid;

#ifdef USE_MPI
	if (mpiMode == MPI_COMM_THREAD) {
		//SCOREP_RECORDING_OFF();
		asagi::Grid::stopCommThread();
		//SCOREP_RECORDING_ON();
	}
#endif // USE_MPI

	logInfo(rank) << "Initializing velocity field. Done.";

#else // USE_ASAGI
	logError() << "This version does not support ASAGI";
#endif // USE_ASAGI
}

}
