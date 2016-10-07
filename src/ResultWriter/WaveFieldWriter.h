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
 */

#ifndef WAVE_FIELD_WRITER_H
#define WAVE_FIELD_WRITER_H

#include "Parallel/MPI.h"

#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

#include "utils/logger.h"

#include "async/Module.h"

#include "Geometry/refinement/VariableSubSampler.h"
#include "Monitoring/instrumentation.fpp"
#include "WaveFieldWriterExecutor.h"

namespace seissol
{

namespace writer
{

class WaveFieldWriter : private async::Module<WaveFieldWriterExecutor, WaveFieldInitParam, WaveFieldParam>
{
	/** True if wave field output is enabled */
	bool m_enabled;

	/** False if entire region is to be written */
	bool m_extractRegion;

	/** The asynchronous executor */
	WaveFieldWriterExecutor m_executor;

	/** Variable buffer ids (high and low order variables) */
	int m_variableBufferIds[2];

	/** The output prefix for the filename */
	std::string m_outputPrefix;

	/** The variable subsampler for the refined mesh */
	refinement::VariableSubsampler<double>* m_variableSubsampler;

	/** Number of variables */
	unsigned int m_numVariables;

	/** Flag indicated which variables should be written */
	bool* m_outputFlags;

	/** Refined number of cells */
	unsigned int m_numCells;

	/** Unrefined (low order) number of cells */
	unsigned int m_numLowCells;

	/** Pointer to the degrees of freedom */
	const double* m_dofs;

	/** Pointer to the plastic strain */
	const double* m_pstrain;

	/** Mapping from the cell order to dofs order */
	unsigned int* m_map;

	/** Time of the last output (makes sure output is not written twice at the end) */
	double m_lastTimeStep;

	/** The tolerance in the time for ignoring duplicate time steps */
	double m_timeTolerance;

	/** The current output time step */
	unsigned int m_timestep;

	/** Checks if a vertex given by the vertexCoords lies inside the boxBounds */
	/*   The boxBounds is in the format: xMin, xMax, yMin, yMax, zMin, zMax */
	bool vertexInBox(const double * const boxBounds, const double * const vertexCoords) {
		if (vertexCoords[0] <= boxBounds[1] && vertexCoords[0] >= boxBounds[0] &&
			vertexCoords[1] <= boxBounds[3] && vertexCoords[1] >= boxBounds[2] &&
			vertexCoords[2] <= boxBounds[5] && vertexCoords[2] >= boxBounds[4]) {
			return true;
		} else {
			return false;
		}
	}

public:
	WaveFieldWriter()
		: m_enabled(false),
		  m_extractRegion(false),
		  m_variableSubsampler(0L),
		  m_numVariables(0),
		  m_outputFlags(0L),
		  m_numCells(0), m_numLowCells(0),
		  m_dofs(0L), m_pstrain(0L),
		  m_map(0L),
		  m_lastTimeStep(-1),
		  m_timeTolerance(0),
		  m_timestep(0)
	{
	}

	/**
	 * Activate the wave field output
	 */
	void enable()
	{
		m_enabled = true;
	}

	/**
	 * @return True if wave field output is enabled, false otherwise
	 */
	bool isEnabled() const
	{
		return m_enabled;
	}

	/**
	 * Set the output prefix for the filename
	 */
	void setFilename(const char* outputPrefix)
	{
		m_outputPrefix = outputPrefix;
	}

	/**
	 * Called by ASYNC on all ranks
	 */
	void setUp()
	{
		setExecutor(m_executor);
	}

	/**
	 * Initialize the wave field ouput
	 *
	 * @param map The mapping from the cell order to dofs order
	 * @param timeTolerance The tolerance in the time for ignoring duplicate time steps
	 */
	void init(unsigned int numVars, int order, int numAlignedDOF,
			const MeshReader &meshReader,
			const double* dofs,  const double* pstrain,
			unsigned int* map,
			int refinement, int timestep, int* outputMask, double* outputRegionBounds,
			double timeTolerance);

	/**
	 * @return The current time step of the wave field output
	 */
	unsigned int timestep() const
	{
		return m_timestep;
	}

	/**
	 * Write a time step
	 */
	void write(double time)
	{
		SCOREP_USER_REGION("WaveFieldWriter_write", SCOREP_USER_REGION_TYPE_FUNCTION);

		if (!m_enabled)
			return;

		const int rank = seissol::MPI::mpi.rank();

		if (time <= m_lastTimeStep + m_timeTolerance) {
			// Ignore duplicate time steps. Might happen at the end of a simulation
			logInfo(rank) << "Ignoring duplicate time step at time " << time;
			return;
		}

		SCOREP_USER_REGION_DEFINE(r_wait);
		SCOREP_USER_REGION_BEGIN(r_wait, "wavfieldwriter_wait", SCOREP_USER_REGION_TYPE_COMMON);
		logInfo(rank) << "Waiting for last wave field.";
		wait();
		SCOREP_USER_REGION_END(r_wait);

		logInfo(rank) << "Writing wave field at time" << utils::nospace <<  time << '.';

		unsigned int nextId = m_variableBufferIds[0];
		for (unsigned int i = 0; i < m_numVariables; i++) {
			if (!m_outputFlags[i])
				continue;

			double* managedBuffer = async::Module<WaveFieldWriterExecutor,
					WaveFieldInitParam, WaveFieldParam>::managedBuffer<double*>(nextId);
			m_variableSubsampler->get(m_dofs, m_map, i, managedBuffer);

			sendBuffer(nextId, m_numCells*sizeof(double));

			nextId++;
		}

		if (m_pstrain) {
			for (unsigned int i = 0; i < WaveFieldWriterExecutor::NUM_LOWVARIABLES; i++) {
				double* managedBuffer = async::Module<WaveFieldWriterExecutor,
						WaveFieldInitParam, WaveFieldParam>::managedBuffer<double*>(m_variableBufferIds[1]+i);

#ifdef _OPENMP
				#pragma omp parallel for schedule(static)
#endif // _OPENMP
				for (unsigned int j = 0; j < m_numLowCells; j++)
					managedBuffer[j] = m_pstrain[m_map[j]
							* WaveFieldWriterExecutor::NUM_LOWVARIABLES + i];

				sendBuffer(m_variableBufferIds[1]+i, m_numLowCells*sizeof(double));
			}
		}

		WaveFieldParam param;
		param.time = time;
		call(param);

		// Update last time step
		m_lastTimeStep = time;
		m_timestep++;

		logInfo(rank) << "Writing wave field at time" << utils::nospace << time << ". Done.";
	}

	/**
	 * Close wave field writer and free resources
	 */
	void close()
	{
		// Cleanup the executor
		if (m_enabled)
			wait();

		finalize();

		if (!m_enabled)
			return;

		delete m_variableSubsampler;
		m_variableSubsampler = 0L;
		delete [] m_outputFlags;
		m_outputFlags = 0L;
		if (m_extractRegion) {
			delete [] m_map;
			m_map = 0L;
		}
	}

	void tearDown()
	{
		if (!m_enabled)
			return;

		m_executor.finalize();
	}
};

}

}

#endif // WAVE_FIELD_WRITER_H
