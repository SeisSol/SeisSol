/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
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

#ifndef ASAGIREADER_H
#define ASAGIREADER_H

#ifdef USE_ASAGI

#include "Parallel/MPI.h"

#include <cassert>
#include <string>

#include <asagi.h>

#include "glm/vec3.hpp"
#include "glm/vector_relational.hpp"
#include "glm/gtc/round.hpp"

#include "utils/env.h"
#include "utils/logger.h"

#include "AsagiModule.h"
#include "Monitoring/instrumentation.fpp"

namespace seissol
{

namespace asagi
{

typedef double vertex_t[3];

enum NUMACache_Mode
{
	NUMA_OFF, NUMA_ON, NUMA_CACHE
};

class AsagiReader
{
private:
	/** Prefix for environment variables */
	const std::string m_envPrefix;

	/** The ASAGI grid */
	::asagi::Grid* m_grid;

	/** Number of threads used by ASAGI */
	unsigned int m_asagiThreads;

	/** Number of values in the grid */
	unsigned int m_numValues;

public:
	AsagiReader(const char* envPrefix)
		: m_envPrefix(envPrefix), m_grid(0L), m_numValues(0)
	{
	}

	~AsagiReader()
	{
		close();
	}

	void open(const char* file, unsigned int numValues)
	{
		assert(numValues > 0);

		m_numValues = numValues;

		SCOREP_USER_REGION("AsagiReader_open", SCOREP_USER_REGION_TYPE_FUNCTION);

		const int rank = seissol::MPI::mpi.rank();

		m_grid = ::asagi::Grid::createArray();

		// Set MPI mode
		if (AsagiModule::mpiMode() != MPI_OFF) {
#ifdef USE_MPI
			m_grid->setComm(seissol::MPI::mpi.comm());
#endif // USE_MPI

			if (AsagiModule::mpiMode() == MPI_COMM_THREAD)
				m_grid->setParam("MPI_COMMUNICATION", "THREAD");
		}

		// Set NUMA mode
		m_asagiThreads = utils::Env::get((m_envPrefix  + "_NUM_THREADS").c_str(), 0u);
		if (m_asagiThreads == 0)
			m_asagiThreads = AsagiModule::totalThreads();
		else if (m_asagiThreads > AsagiModule::totalThreads()) {
			logWarning(rank) << "Only" << AsagiModule::totalThreads()
					<< "threads can be used for ASAGI initialization.";
			m_asagiThreads = AsagiModule::totalThreads();
		}

		if (AsagiModule::mpiMode() == MPI_COMM_THREAD)
			m_asagiThreads--; // one thread is used for communication

		m_grid->setThreads(m_asagiThreads);

		switch (getNUMAMode()) {
		case NUMA_ON:
			m_grid->setParam("NUMA_COMMUNICATION", "ON");
			break;
		case NUMA_OFF:
			m_grid->setParam("NUMA_COMMUNICATION", "OFF");
			break;
		case NUMA_CACHE:
			m_grid->setParam("NUMA_COMMUNICATION", "CACHE");
			break;
		}

		// Set vertex centered grid
		m_grid->setParam("VALUE_POSITION", "VERTEX_CENTERED");

		// Set additional parameters
		std::string blockSize = utils::Env::get((m_envPrefix  + "_BLOCK_SIZE").c_str(), "64");
		m_grid->setParam("BLOCK_SIZE_0", blockSize.c_str());
		m_grid->setParam("BLOCK_SIZE_1", blockSize.c_str());
		m_grid->setParam("BLOCK_SIZE_2", blockSize.c_str());

		std::string cacheSize = utils::Env::get((m_envPrefix  + "_CACHE_SIZE").c_str(), "128");
		m_grid->setParam("CACHE_SIZE", cacheSize.c_str());

		m_grid->setParam("VARIABLE", "data");

		// Read the data
		//SCOREP_RECORDING_OFF();
#ifdef _OPENMP
		#pragma omp parallel num_threads(m_asagiThreads)
#endif // _OPENMP
		{
			::asagi::Grid::Error err = m_grid->open(file);
			if (err != ::asagi::Grid::SUCCESS)
				logError() << "Could not open ASAGI grid:" << err;
		}
		//SCOREP_RECORDING_ON();

		if (m_grid->getVarSize() != numValues * sizeof(float))
			logError() << "Invalid variable size in material file";
	}

	/**
	 * @param[out] outside Number of cells found outside the ASAGI grid
	 */
	template<class Setter>
	void read(const vertex_t* coords, unsigned int numElements,
			Setter setter, const float* defaultValues,
			unsigned long &outside)
	{
		SCOREP_USER_REGION("AsagiReader_read", SCOREP_USER_REGION_TYPE_FUNCTION);

		const int rank = seissol::MPI::mpi.rank();

		// Grid dimensions
		const glm::dvec3 min(m_grid->getMin(0), m_grid->getMin(1), m_grid->getMin(2));
		const glm::dvec3 max(m_grid->getMax(0), m_grid->getMax(1), m_grid->getMax(2));
		const glm::dvec3 delta(m_grid->getDelta(0), m_grid->getDelta(1), m_grid->getDelta(2));

		unsigned long outsideInt = 0;

		// Allocate buffer
		float* values[8];
		for (unsigned int i = 0; i < 8; i++)
			values[i] = new float[m_numValues];
		float* interpolValues = new float[m_numValues];

		//SCOREP_RECORDING_OFF();
#ifdef _OPENMP
		#pragma omp parallel for num_threads(m_asagiThreads) reduction(+: outsideInt)
#endif // _OPENMP
		for (int i = 0; i < numElements; i++) {
			const glm::dvec3 coord(coords[i][0], coords[i][1], coords[i][2]);

			// Compute the surrounding coordinates
			glm::dvec3 shiftedBary = coord - min;
			glm::dvec3 lowCoord = glm::floorMultiple(shiftedBary, delta) + min;
			glm::dvec3 highCoord = glm::ceilMultiple(shiftedBary, delta) + min;

			// Fix low/high if they are outside of the domain (but the bary center is inside)
			// -> Should not be necessary if we use vertex centered grid
//			if (glm::all(glm::greaterThanEqual(coord, min))
//				&& glm::any(glm::lessThan(lowCoord, min)))
//				lowCoord = glm::max(lowCoord, min);
//			if (glm::all(glm::lessThanEqual(coord, max))
//				&& glm::any(glm::greaterThan(highCoord, max)))
//				highCoord = glm::min(highCoord, max);

			if (glm::any(glm::lessThan(lowCoord, min))
				|| glm::any(glm::greaterThan(highCoord, max))) {
				// Outside the ASAGI domain

				setter.set(i, defaultValues);

				outside++;
			} else {
				// Get all 8 values
				double pos[3];
				pos[0] = lowCoord.x; pos[1] = lowCoord.y; pos[2] = lowCoord.z;
				m_grid->getBuf(values[0], pos);
				pos[0] = lowCoord.x; pos[1] = lowCoord.y; pos[2] = highCoord.z;
				m_grid->getBuf(values[1], pos);
				pos[0] = lowCoord.x; pos[1] = highCoord.y; pos[2] = lowCoord.z;
				m_grid->getBuf(values[2], pos);
				pos[0] = lowCoord.x; pos[1] = highCoord.y; pos[2] = highCoord.z;
				m_grid->getBuf(values[3], pos);
				pos[0] = highCoord.x; pos[1] = lowCoord.y; pos[2] = lowCoord.z;
				m_grid->getBuf(values[4], pos);
				pos[0] = highCoord.x; pos[1] = lowCoord.y; pos[2] = highCoord.z;
				m_grid->getBuf(values[5], pos);
				pos[0] = highCoord.x; pos[1] = highCoord.y; pos[2] = lowCoord.z;
				m_grid->getBuf(values[6], pos);
				pos[0] = highCoord.x; pos[1] = highCoord.y; pos[2] = highCoord.z;
				m_grid->getBuf(values[7], pos);

				// Do trilinear interpolation:
				// https://en.wikipedia.org/wiki/Trilinear_interpolation
				glm::dvec3 d = (coord - lowCoord) / delta;

				for (unsigned int j = 0; j < 3; j++) {
					double c0[4];
					for (unsigned int k = 0; k < 4; k++)
						c0[k] = values[k][j]*(1-d.x) + values[k+4][j]*d.x;

					double c1[2];
					for (unsigned int k = 0; k < 2; k++)
						c1[k] = c0[k]*(1-d.y) + c0[k+1]*d.y;

					interpolValues[j] = c1[0]*(1-d.z) + c1[1]*d.z;
				}

				setter.set(i, interpolValues);
			}
		}
		//SCOREP_RECORDING_ON();

		for (unsigned int i = 0; i < 8; i++)
			delete [] values[i];
		delete [] interpolValues;

#ifdef USE_MPI
		if (rank == 0)
			MPI_Reduce(MPI_IN_PLACE, &outsideInt, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, seissol::MPI::mpi.comm());
		else
			MPI_Reduce(&outsideInt, 0L, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, seissol::MPI::mpi.comm());
#endif // USE_MPI
		outside = outsideInt;
	}

	void close()
	{
		delete m_grid;
		m_grid = 0L;
	}

private:
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
};

}

}

#endif // USE_ASAGI

#endif // ASAGIREADER_H